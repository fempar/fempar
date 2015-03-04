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
  use types
  use memor
  use fem_matrix_names
  use fem_partition_names
  use array_names
  use stdio
  use psb_penv_mod
#ifdef memcheck
  use iso_c_binding
#endif
  
  ! Module associated with the F90 interface to Trilinos.
  ! Remember: the F90 interface to Trilinos requires C
  ! interoperability (i.e., iso_c_binding module)
  !use for_trilinos_shadow_interfaces

  ! Parallel modules
  use par_graph_names
  use par_context_names
  use par_partition_names

  implicit none
# include "debug.i90"

  private

  type par_matrix
     !type( epetra_crsmatrix ) :: epm

     ! Data structure which stores the local part 
     ! of the matrix mapped to the current processor.
     ! This is required for both eb and vb data 
     ! distributions
     type( fem_matrix )       :: f_matrix

     type(par_graph), pointer :: &
        p_graph => NULL()           ! Associated par_graph
     
     type(par_partition), pointer :: &
        p_part => NULL()            ! Associated (ROW) par_partition
     
     type(par_partition), pointer :: &
        p_part_cols => NULL()            ! Associated (COL) par_partition

  end type par_matrix

  interface par_matrix_free
     module procedure par_matrix_free_one_shot, par_matrix_free_progressively
  end interface par_matrix_free

  ! Types
  public :: par_matrix

  ! Functions
  public :: par_matrix_create, par_matrix_graph, par_matrix_fill_val, &
         &  par_matrix_alloc, par_matrix_free,    & 
         &  par_matrix_assembly, par_matrix_info, &
         &  par_matrix_print, par_matrix_print_matrix_market, &
         &  par_matrix_zero, &
         &  par_matrix_bcast, par_matrix_fine_task

!***********************************************************************
! Allocatable arrays of type(par_matrix)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(par_matrix)
# define var_size 80
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine par_matrix_create(storage,type,symm,nd1,nd2,p_part,p_part_cols,p_matrix,def)
    implicit none
    integer(ip)     , intent(in)            :: storage,type,symm,nd1,nd2
    type(par_partition), target, intent(in) :: p_part
    type(par_partition), target, intent(in) :: p_part_cols
    type(par_matrix), intent(out)           :: p_matrix
    integer(ip)     , optional, intent(in)  :: def

    call fem_matrix_create(storage,type,symm,nd1,nd2,p_matrix%f_matrix,def)
    p_matrix%p_part      => p_part 
    p_matrix%p_part_cols => p_part_cols 

  end subroutine par_matrix_create

  subroutine par_matrix_graph(p_graph,p_matrix)
    implicit none
    type(par_graph) , target, intent(in)   :: p_graph
    type(par_matrix), intent(inout)        :: p_matrix

    ! Pointer to part/context object is required
    assert ( associated(p_graph%p_part) )
    assert ( associated(p_graph%p_part%p_context) )
    assert ( p_graph%p_part%p_context%created .eqv. .true.)
    
    ! Trilinos requires (unsymm.) CSR sparse matrix format
    assert ( (.not. (p_graph%p_part%p_context%handler == trilinos)) .or. (p_matrix%f_matrix%type == csr_mat .and. p_matrix%f_matrix%symm==symm_false) )

    ! Point to input target parallel graph
    p_matrix%p_graph => p_graph

    if(p_graph%p_part%p_context%iam>=0) call fem_matrix_graph ( p_graph%f_graph, p_matrix%f_matrix )

  end subroutine par_matrix_graph


  subroutine par_matrix_fill_val(p_matrix)
    implicit none
    type(par_matrix), intent(inout)          :: p_matrix

    if(p_matrix%p_part%p_context%iam>=0) then
       call fem_matrix_fill_val(p_matrix%f_matrix)
    else
       return
    end if

    ! write (*,*) size(p_matrix%f_matrix%a,1), size(p_matrix%f_matrix%a,2), size(p_matrix%f_matrix%a,3) DBG:

    if ( p_matrix%p_part%p_context%handler == trilinos ) then 
       ! epetra_crsmatrix (i.e., scalar matrix)
       if (p_matrix%f_matrix%storage == scal ) then 
         call par_matrix_epetra_crsmatrix_create (p_matrix)
       else ! epetra_vbrmatrix (i.e., block matrix)
         write (0,*) 'Error: trilinos shadow interfaces do not yet support epetra_vbrmatrix object (i.e., block matrices)'
         stop
       end if 
    end if

  end subroutine par_matrix_fill_val

  subroutine par_matrix_alloc(storage,type,symm,nd1,nd2,p_graph,p_matrix,def)
    implicit none

    integer(ip)     , intent(in)           :: storage,type,symm,nd1,nd2
    type(par_graph) , target, intent(in)   :: p_graph
    type(par_matrix), intent(out)          :: p_matrix
    integer(ip)     , optional, intent(in) :: def

    call par_matrix_create(storage,type,symm,nd1,nd2,p_graph%p_part,p_graph%p_part_cols,p_matrix,def)
    call par_matrix_graph(p_graph,p_matrix)
    call par_matrix_fill_val(p_matrix)

    ! ! Pointer to part/context object is required
    ! assert ( associated(p_graph%p_part) )
    ! assert ( associated(p_graph%p_part%p_context) )
    ! assert ( p_graph%p_part%p_context%created .eqv. .true.)
    
    ! ! Trilinos requires (unsymm.) CSR sparse matrix format
    ! assert ( (.not. (p_graph%p_part%p_context%handler == trilinos)) .or. (type == csr_mat .and. symm==symm_false) )

    ! ! Point to input target parallel graph
    ! p_matrix%p_graph => p_graph


    ! if(p_graph%p_part%p_context%iam>=0) then
    !    ! Allocate local part
    !    if(present(def)) then
    !       !write(*,*) 'Present',def
    !       call fem_matrix_alloc ( storage, type, symm, nd1, nd2, & 
    !            &   p_graph%f_graph, p_matrix%f_matrix, def )
    !    else
    !       !write(*,*) 'not present'
    !       call fem_matrix_alloc ( storage, type, symm, nd1, nd2, & 
    !            &   p_graph%f_graph, p_matrix%f_matrix )
    !    end if
    ! else
    !    p_matrix%f_matrix%nd1     = nd1
    !    p_matrix%f_matrix%nd2     = nd2
    !    p_matrix%f_matrix%symm    = symm
    !    if (present(def)) then
    !       p_matrix%f_matrix%sign = def
    !    else
    !       p_matrix%f_matrix%sign = unknown
    !    end if
    !    p_matrix%f_matrix%storage = storage 
    !    p_matrix%f_matrix%type    = type
    ! end if

    ! if(p_graph%p_part%p_context%iam<0) return

    ! ! write (*,*) size(p_matrix%f_matrix%a,1), size(p_matrix%f_matrix%a,2), size(p_matrix%f_matrix%a,3) DBG:

    ! if ( p_graph%p_part%p_context%handler == trilinos ) then 
    !    ! epetra_crsmatrix (i.e., scalar matrix)
    !    if ( storage == scal ) then 
    !      call par_matrix_epetra_crsmatrix_create (p_matrix)
    !    else ! epetra_vbrmatrix (i.e., block matrix)
    !      write (0,*) 'Error: trilinos shadow interfaces do not yet support epetra_vbrmatrix object (i.e., block matrices)'
    !      stop
    !    end if 
    ! end if

  end subroutine par_matrix_alloc

  !=============================================================================
  subroutine par_matrix_free_one_shot(p_matrix)
    implicit none

    type(par_matrix), intent(inout) :: p_matrix
    call par_matrix_free_progressively(p_matrix, free_only_values)
    call par_matrix_free_progressively(p_matrix, free_only_struct)
    call par_matrix_free_progressively(p_matrix, free_clean)
  end subroutine par_matrix_free_one_shot

  !=============================================================================
  subroutine par_matrix_free_progressively(p_matrix, mode)
    implicit none
    type(par_matrix), intent(inout) :: p_matrix
    integer(ip)     , intent(in)    :: mode

    ! The routine requires the partition/context info
    assert ( associated(p_matrix%p_part) )
    assert ( associated(p_matrix%p_part%p_context) )
    assert ( p_matrix%p_part%p_context%created .eqv. .true.)
    assert ( mode == free_clean .or. mode == free_only_struct .or. mode == free_only_values )

    if(p_matrix%p_part%p_context%iam<0) return

    if ( mode == free_clean ) then
       nullify ( p_matrix%p_part )
       nullify ( p_matrix%p_part_cols )
    else if ( mode == free_only_struct ) then
       ! AFM: This nullification cannot be here as this means that it will not be longer possible
       !      to access p_matrix%p_part after "free_only_struct"ing a par_matrix
       !      (and it is done in many parts of the code). I will move it to free_clean.
       ! AFM: The comment above NO longer applies as par_matrix now directly points to
       !      the par_partition instance 
       nullify ( p_matrix%p_graph )
       ! else if ( mode == free_only_values ) then
       !    if ( p_matrix%p_part%p_context%handler == trilinos  ) then 
       !       ! epetra_crsmatrix (i.e., scalar matrix)
       !       if ( p_matrix%f_matrix%storage == scal ) then 
       !          call epetra_crsmatrix_destruct ( p_matrix%epm )
       !       else ! epetra_vbrmatrix (i.e., block matrix)
       !          write (0,*) 'Error: trilinos shadow interfaces do not yet support epetra_vbrmatrix object (i.e., block matrices)'
       !          stop
       !       end if
       !    end if
    end if

    ! Free local part
    call fem_matrix_free ( p_matrix%f_matrix, mode )

  end subroutine par_matrix_free_progressively
  !=============================================================================
  subroutine par_matrix_assembly(nn,ln,ea,p_mat)
    implicit none
    
    ! Parameters 
    integer(ip) , intent(in)        :: nn
    integer(ip) , intent(in)        :: ln(:)
    type(array_rp2) , intent(in)        :: ea
    type(par_matrix), intent(inout) :: p_mat
     
    ! TO-DO: currently we are calling serial matrix assembly
    ! routine for both element-based and vertex-based
    ! partitionings. For vertex-based partitionings this
    ! is not very intelligent, as it involves assembling
    ! (with partial contributions) the rows corresponding 
    ! to external nodes. Rows corresponding to external nodes
    ! are just wasted by the underlying linear algebra package.
    ! It would be more intelligent to reduce memory traffic to
    ! avoid accessing (i.e., read+write) these rows by splitting
    ! the external loop on the elements into two parts (internal+boundary
    ! elements)
    call fem_matrix_assembly ( nn, ln, ea, p_mat%f_matrix )

  end subroutine par_matrix_assembly


  !===================================================
  ! Creates and fills epetra_crsmatrix p_matrix%epm
  ! as a view of the local matrix f_matrix
  !===================================================
  ! subroutine par_matrix_epetra_crsmatrix_create (p_matrix)
!     implicit none

!     ! Parameters 
!     type(par_matrix), intent(inout) :: p_matrix

!     ! Local variables
!     integer (c_int)                 :: ierrc

!     ! epetra_crsgraph is required
!     assert ( associated(p_matrix%p_graph) )

!     ! Vertex-based partitioning/trilinos handler required for Epetra
!     assert ( p_matrix%p_part%f_part%ptype == vertex_based )
!     assert ( p_matrix%p_part%p_context%handler == trilinos )

!     call epetra_crsmatrix_construct ( p_matrix%epm, p_matrix%p_graph%epg )

!     ! Insert row/column entries information via calls to InsertMyValues
!     call insert_my_values ( p_matrix%epm,                                                   &  
!          (p_matrix%p_part%f_part%nmap%ni +                       & 
!          p_matrix%p_part%f_part%nmap%nb)*p_matrix%f_matrix%nd1,  &
!          p_matrix%p_graph%f_graph%ia,                                    &
!          p_matrix%p_graph%f_graph%ja,                                    &
!          p_matrix%f_matrix%a )

!     ! Signal to inform that we are done with matrix entries insertion
!     ! call epetra_crsmatrix_fill_complete ( p_matrix%epm, ierrc )

!     ! Signal to inform that we are done                       
! !!$    call epetra_crsmatrix_fill_complete ( p_matrix%epm, p_matrix%p_part%row_map(p_matrix%f_matrix%nd2), &
! !!$         & p_matrix%p_part%row_map(p_matrix%f_matrix%nd1), ierrc );
! !!$    assert (ierrc == 0)

!   end subroutine par_matrix_epetra_crsmatrix_create

  ! Private routine which unpacks allocatable derived data type 
  ! members as explicit size arrays
  ! subroutine insert_my_values (epm, n, ia, ja, a)
  !   ! Parameters
  !   type(epetra_crsmatrix)  , intent(inout) :: epm
  !   integer (c_int)         , intent(in)    :: n
  !   integer (c_int)         , intent(in)    :: ia(n+1)
  !   integer (c_int)         , intent(in)    :: ja(ia(n+1)-1)
  !   real (c_double)         , intent(in)    :: a (ia(n+1)-1)

  !   ! Local variables
  !   integer (c_int)                 :: i, ierrc

  !   ! Insert row and column entries information via calls to InsertMyValues
  !   do i=1, n
  !       call epetra_crsmatrix_insert_my_values ( epm, i-1,              & 
  !            &                                   ia(i+1)-ia(i),         & 
  !            &                                   a(ia(i):(ia(i+1)-1)),  & 
  !            &                                   ja(ia(i):(ia(i+1)-1)), & 
  !            &                                   ierrc )
  !       assert (ierrc == 0)
  !   end do 
  
  ! end subroutine insert_my_values

  subroutine par_matrix_info ( p_mat, me, np )
    implicit none

    ! Parameters 
    type(par_matrix), intent(in)    :: p_mat
    integer         , intent(out)   :: me
    integer         , intent(out)   :: np

   ! Pointer to part/context object is required
    assert ( associated(p_mat%p_part) )
    assert ( associated(p_mat%p_part%p_context) )
    assert ( p_mat%p_part%p_context%created .eqv. .true.)
    me = p_mat%p_part%p_context%iam
    np = p_mat%p_part%p_context%np

    !call par_context_info ( p_mat%p_part%p_context, me, np )

  end subroutine par_matrix_info

  !=============================================================================
  subroutine par_matrix_print(lunou, p_matrix)
    implicit none
    type(par_matrix)  ,  intent(in) :: p_matrix
    integer(ip)      ,  intent(in) :: lunou

    ! p_graph%p_part is required within this subroutine
    assert ( associated(p_matrix%p_part) )
    
    ! p_graph%p_part%p_context is required within this subroutine
    assert ( associated(p_matrix%p_part%p_context) )
    
    ! if ( p_matrix%p_part%p_context%handler == trilinos ) then 
    !     call epetra_crsmatrix_print ( p_matrix%epm )
    ! end if

  end subroutine par_matrix_print

  subroutine par_matrix_print_matrix_market ( dir_path, prefix, p_mat, global )
    implicit none
    ! Parameters
    character *(*)  , intent(in)           :: dir_path
    character *(*)  , intent(in)           :: prefix
    type(par_matrix), intent(in)           :: p_mat
    logical         , intent(in), optional :: global  

    ! Locals
    integer         :: iam, num_procs, lunou
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id
    logical         :: global_

    assert ( associated(p_mat%p_part%p_context) )
    assert ( p_mat%p_part%p_context%created .eqv. .true.)
    if(p_mat%p_part%p_context%iam<0) return

    if (present(global)) then
      global_ = global
    else
      global_ = .false.
    end if

    name = trim(prefix) // '.par_matrix' // '.mtx'

    ! Get context info
    call par_context_info ( p_mat%p_part%p_context, iam, num_procs )

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

    ! Read fem_partition data from path_file file
    lunou =  io_open (trim(dir_path) // '/' // trim(name) // '.' // trim(zeros) // trim(part_id), 'write')

    if (global_) then  
       call fem_matrix_print_matrix_market ( lunou, p_mat%f_matrix,               &
                                             p_mat%p_part%f_part%nmap%ng, &
                                             p_mat%p_part%f_part%nmap%l2g )
    else
       call fem_matrix_print_matrix_market ( lunou, p_mat%f_matrix )
    end if 

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
    type(par_matrix), intent(inout)    :: p_matrix

    ! p_part%p_context is required within this subroutine
    assert ( associated(p_matrix%p_part%p_context) )
    assert ( p_matrix%p_part%p_context%created .eqv. .true.)

    if(p_matrix%p_part%p_context%iam<0) return

    call fem_matrix_zero ( p_matrix%f_matrix )
  end subroutine par_matrix_zero

  subroutine par_matrix_bcast (p_matrix, conv)
    implicit none
    ! Parameters 
    type(par_matrix), target, intent(in)    :: p_matrix    
    logical                 , intent(inout) :: conv
    assert ( associated(p_matrix%p_part   ) )
    call par_partition_bcast(p_matrix%p_part,conv)
  end subroutine par_matrix_bcast

  function par_matrix_fine_task (p_matrix)
    implicit none
    logical                         :: par_matrix_fine_task
    type(par_matrix), intent(in)    :: p_matrix

    assert ( associated(p_matrix%p_part   ) )
    assert ( associated(p_matrix%p_part%p_context) ) 
    assert ( p_matrix%p_part%p_context%created .eqv. .true.)

    par_matrix_fine_task = .true. 

    ! Only if global context has been created (i.e.,
    ! there are coarse MPI Tasks) it is required to
    ! broadcast conv
    if ( associated(p_matrix%p_part%g_context) ) then
       assert ( p_matrix%p_part%g_context%created .eqv. .true.)
       par_matrix_fine_task = (p_matrix%p_part%p_context%iam >= 0)
    end if

  end function par_matrix_fine_task


end module par_matrix_names
