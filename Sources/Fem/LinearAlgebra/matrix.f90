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
module matrix_names
  use types_names
  use memor_names
  use sort_names
  use graph_names
  use vector_names
  use matvec_dof_names

  ! Abstract types
  use base_operand_names
  use base_operator_names

#ifdef memcheck
use iso_c_binding
#endif

  implicit none
  private

  ! Constants:
  ! a) matrix_t types
  integer(ip), parameter :: csr_mat=10
  ! b) matrix symmetry
  integer(ip), parameter :: symm_true =0
  integer(ip), parameter :: symm_false=1
  ! c) matrix sign
  integer(ip), parameter :: positive_definite     = 0
  integer(ip), parameter :: positive_semidefinite = 1
  integer(ip), parameter :: indefinite            = 2 ! Both positive and negative eigenvalues
  integer(ip), parameter :: unknown               = 3 ! No info

  ! Matrix
  type, extends(base_operator_t) :: matrix_t
     integer(ip)                :: &
!!!          storage=undef_sto,         &      ! Storage layout (blk: block; scal: scalar)
          symm=symm_false,           &         ! Flag for symmetry
          sign=positive_definite,    &         ! Flag for positiveness
          type=csr_mat                         ! matrix_t type


     ! We need to decide how blocks will be stored (transposed or not)
     ! Currently in felap u is transposed but l is not. Seems better to
     ! transpose both. What about trilinos?
     real(rp), allocatable      :: &
          !d(:),                  &         ! Diagonal components (neq) 
          !l(:),                  &         ! Lower in css (nzs)        
          !u(:),                  &         ! Upper in css (nzs)        
          a(:)                             ! Lower/Upper part components ordered as:
     ! Rows (nzt)                
     ! Cols (nzt)                

     type(graph_t),pointer :: &
          gr => NULL()                      ! Associated graph

   contains
     procedure  :: apply     => matrix_apply
     procedure  :: apply_fun => matrix_apply_fun
     procedure  :: free      => matrix_free_tbp
  end type matrix_t

  interface matrix_free
     module procedure matrix_free_one_shot, matrix_free_progressively
  end interface matrix_free

  ! Constants
  public :: csr_mat, symm_true, symm_false
  public :: positive_definite, positive_semidefinite, indefinite, unknown

  ! Types
  public :: matrix_t

  ! Functions
  public :: matrix_create, matrix_graph, matrix_fill_val, &
       &    matrix_alloc, matrix_copy, matrix_print, matrix_read,    &
       &    matrix_print_matrix_market,  & 
       &    matrix_free, matrix_info, matrix_zero,  &
       &    matrix_sum, matrix_transpose, &
       &    matrix_compose_name_matrix_market, matrix_read_matrix_market, &
       &    matvec, matvec_trans, matmat, matmat_trans

!***********************************************************************
! Allocatable arrays of type(matrix_t)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(matrix_t)
# define var_size 80
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine matrix_create(type,symm,mat,def)
    implicit none
    integer(ip)     , intent(in)           :: type, symm
    type(matrix_t), intent(out)          :: mat
    integer(ip)     , optional, intent(in) :: def

    assert ( type == csr_mat )
    assert ( symm == symm_true .or. symm == symm_false )

    mat%symm    =  symm    ! Flag for symmetry
    mat%type    =  type    ! matrix_t type (csr_mat)

    mat%sign = unknown
    if(present(def)) then
       assert ( def == positive_definite .or. def == positive_semidefinite  .or. def == indefinite .or. def == unknown )
       mat%sign = def
    end if

  end subroutine matrix_create

  subroutine matrix_graph(gr,mat)
    implicit none
    type(graph_t) , target, intent(in) :: gr
    type(matrix_t), intent(inout)      :: mat
    mat%gr => gr
  end subroutine matrix_graph

  subroutine matrix_fill_val(mat)
    implicit none
    type(matrix_t), intent(inout)      :: mat
    integer(ip)                          :: neq, nzs, nzt

    neq  =  mat%gr%nv              ! Number of rows and columns (equations)

    nzt = mat%gr%ia(mat%gr%nv+1)-1
    call memalloc(nzt,mat%a,__FILE__,__LINE__)
    mat%a = 0.0_rp

  end subroutine matrix_fill_val

  subroutine matrix_alloc(type,symm,gr,mat,def)
    implicit none
    integer(ip)     , intent(in)           :: type, symm
    type(graph_t) , target, intent(in)   :: gr
    type(matrix_t), intent(out)          :: mat
    integer(ip)     , optional, intent(in) :: def

    call matrix_create(type,symm,mat,def)
    call matrix_graph(gr,mat)
    call matrix_fill_val(mat)

  end subroutine matrix_alloc

  subroutine matrix_copy (imatrix, omatrix)
    implicit none
    ! Parameters 
    type(matrix_t), intent(in)    :: imatrix
    type(matrix_t), intent(inout) :: omatrix

    ! *** IMPORTANT NOTE: This routine assumes that omatrix 
    ! already has an associated graph in omatrix%gr

    ! This routine is only provided for matrices stored in CSR
    ! format. If you need to copy a matrix in a different format,
    ! please extend this routine
    assert ( imatrix%type == csr_mat )

    call matrix_alloc ( imatrix%type, imatrix%symm, omatrix%gr, omatrix, imatrix%sign)
    omatrix%a = imatrix%a

  end subroutine matrix_copy

  !=============================================================================
  subroutine matrix_print(lunou, f_matrix)
    implicit none
    type(matrix_t), intent(in)    :: f_matrix
    integer(ip)     , intent(in)    :: lunou
    integer(ip)                     :: i

    assert ( associated(f_matrix%gr) )
    call graph_print (lunou, f_matrix%gr)

    write (lunou, '(a)')     '*** begin matrix data structure ***'
    if ( f_matrix%type == csr_mat ) then
       ! write (lunou, '(4E20.13)')  f_matrix%a(:, :, :)
       do i=1,f_matrix%gr%nv
          !write(lunou,'(10(e14.7,1x))') f_matrix%a(:,:,f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
          write(lunou,'(10(e25.16,1x))') f_matrix%a(f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
       end do
    end if

  end subroutine matrix_print

  !=============================================================================
  subroutine matrix_read(lunin, f_matrix)
    implicit none
    ! Parameters
    integer(ip)     , intent(in)    :: lunin
    type(matrix_t), intent(inout) :: f_matrix
    
    ! Locals
    integer(ip)                     :: i,nzt

    ! AFM: the following sentence is now NO longer permitted.
    !      f_matrix should be passed to matrix_read in
    !      such a way that the graph is already associated 
    ! allocate(f_matrix%gr)
    ! call graph_read (lunin, f_matrix%gr)
    assert ( associated(f_matrix%gr) )

    if ( f_matrix%type == csr_mat ) then
       nzt = f_matrix%gr%ia(f_matrix%gr%nv+1)-1
       call memalloc(nzt,f_matrix%a,__FILE__,__LINE__)

       do i=1,f_matrix%gr%nv/10
          read(lunin,'(10e15.7)') f_matrix%a(f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
       end do
       if(mod(f_matrix%gr%nv,10)>0) &
            &read(lunin,'(10e15.7)') f_matrix%a(f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
    end if

  end subroutine matrix_read

  !=============================================================================
  subroutine matrix_compose_name_matrix_market ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.mtx'
  end subroutine matrix_compose_name_matrix_market

  subroutine matrix_print_matrix_market (lunou, f_matrix, ng, l2g)
    implicit none
    type(matrix_t), intent(in)           :: f_matrix
    integer(ip)     , intent(in)           :: lunou
    integer(ip)     , intent(in), optional :: ng
    integer(ip)     , intent(in), optional :: l2g (*)

    integer(ip) :: i, j
    integer(ip) :: nv1_, nv2_, nl_


    if ( present(ng) ) then 
       nv1_ = ng
       nv2_ = ng
    else
       nv1_ = f_matrix%gr%nv
       nv2_ = f_matrix%gr%nv2
    end if

    write (lunou,'(a)') '%%MatrixMarket matrix coordinate real general'
    if (f_matrix%gr%type == csr) then
       write (lunou,*) nv1_,nv2_,f_matrix%gr%ia(f_matrix%gr%nv+1)-1
       do i=1,f_matrix%gr%nv
          do j=f_matrix%gr%ia(i),f_matrix%gr%ia(i+1)-1
             if (present(l2g)) then
                write(lunou,'(i12, i12, e32.25)') l2g(i), l2g(f_matrix%gr%ja(j)), f_matrix%a(j)
             else
                write(lunou,'(i12, i12, e32.25)') i, f_matrix%gr%ja(j), f_matrix%a(j)
             end if
          end do
       end do
    else if (f_matrix%gr%type == csr_symm) then
       write (lunou,*) nv1_,nv2_,& 
            2*(f_matrix%gr%ia(f_matrix%gr%nv+1)-1) - f_matrix%gr%nv

       do i=1,f_matrix%gr%nv
          do j=f_matrix%gr%ia(i),f_matrix%gr%ia(i+1)-1
!!$                   if ( j == f_matrix%gr%ia(i) ) then
!!$                      if ( i /= f_matrix%gr%ja(j) ) write (*,*) 'ERR', i, f_matrix%gr%ja(f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1), f_matrix%gr%ia(i)
!!$                      assert ( i ==  f_matrix%gr%ja(j) )
!!$                   end if

             if (present(l2g)) then
                write(lunou,'(i12, i12, e32.25)') l2g(i), l2g(f_matrix%gr%ja(j)), f_matrix%a(j)
             else
                write(lunou,'(i12, i12, e32.25)') i, f_matrix%gr%ja(j), f_matrix%a(j)
             end if
             if (i /= f_matrix%gr%ja(j)) then
                if (present(l2g)) then
                   write(lunou,'(i12, i12, e32.25)') l2g(f_matrix%gr%ja(j)), l2g(i), f_matrix%a(j)
                else
                   write(lunou,'(i12, i12, e32.25)') f_matrix%gr%ja(j), i, f_matrix%a(j)
                end if
             end if
          end do
       end do
    end if

  end subroutine matrix_print_matrix_market

  subroutine matrix_read_matrix_market (lunou, mat, gr, symm, sign)
    implicit none
    integer(ip)     , intent(in)            :: lunou
    type(matrix_t), intent(inout)         :: mat
    type(graph_t) , target, intent(inout) :: gr 
    integer(ip)     , intent(in), optional  :: symm
    integer(ip)     , intent(in), optional  :: sign

    integer(ip) :: i, j, iw1(2),iw2(2)
    integer(ip) :: nv, nw, nz

    integer(ip), allocatable :: ija_work(:,:), ija_index(:)
    real(rp)   , allocatable :: a_work(:)

    mat%type    =  csr_mat ! matrix_t type (csr_mat)

    ! AFM: the following sentence is now NO longer permitted.
    !      Now both mat and graph are passed, and this subroutine is
    !      responsible for both reading the graph and the matrix. 
    !      On return, mat%gr will point to the graph
    ! allocate(mat%gr)

    mat%symm    =  symm_false    ! Flag for symmetry
    if(present(symm)) then
       assert ( symm == symm_true .or. symm == symm_false )
       mat%symm = symm
    end if

    mat%sign = unknown
    if(present(sign)) then
       assert ( sign == positive_definite .or. sign == positive_semidefinite  .or. sign == indefinite .or. sign == unknown )
       mat%sign = sign
    end if

    read (lunou,*) 
    read (lunou,*) nv,nw,nz

    if(nv/=nw) then 
       write (0,*) 'Error: only square matrices are allowed in fempar when reading from matrix market'
       stop
    end if

    call memalloc ( 2, nz, ija_work, __FILE__,__LINE__ )
    call memalloc (    nz, ija_index, __FILE__,__LINE__ )
    call memalloc (    nz, a_work,__FILE__,__LINE__)

    if(mat%symm==symm_false) then
       gr%type = csr
       do i=1,nz
          ! read(lunou,*, end=10,err=10) ija_work(1,i),ija_work(2,i),a_work(i)
          read(lunou,'(i12, i12, e32.25)',end=10, err=10) ija_work(1,i), ija_work(2,i), a_work(i)
          ija_index(i)=i
       end do
    else
       gr%type = csr_symm
       j=1
       do i=1,nz
          read(lunou,'(i12, i12, e32.25)', end=10,err=10) ija_work(1,j),ija_work(2,j),a_work(j)
          ija_index(j)=j
          if(ija_work(1,j)<=ija_work(2,j)) j=j+1 ! Keep it (but if it is the last keep j)
       end do
       nz=j-1 !-nv
       write(*,*) nz
    end if

    ! Sort by ia and ja
    call intsort(2,2,nz,ija_work,ija_index,iw1,iw2)

    ! Allocate graph 
    gr%nv   = nv
    gr%nv2  = nv
    gr%nzt  = nz
    call memalloc ( nv+1 , gr%ia, __FILE__,__LINE__ )
    call memalloc ( nz , gr%ja, __FILE__,__LINE__ )

    ! Count columns on each row
    gr%ia = 0
    do i=1,nz
       gr%ia(ija_work(1,i)+1) = gr%ia(ija_work(1,i)+1) + 1
    end do
    ! Compress
    gr%ia(1)=1
    do i=1,nv
       gr%ia(i+1) = gr%ia(i+1) +gr%ia(i)
    end do
    write(*,*) gr%ia(nv+1)
    ! Copy ja
    gr%ja = ija_work(2,:)
    call memfree(ija_work,__FILE__,__LINE__)

    ! Reorder a
    call memalloc (nz,mat%a,__FILE__,__LINE__)
    do i=1,nz
       mat%a(i)=a_work(ija_index(i))
    end do
    call memfree ( ija_index,__FILE__,__LINE__)
    call memfree (    a_work,__FILE__,__LINE__)

    mat%gr => gr

    return

10  write (0,*) 'Error reading matrix eof or err'
    stop

  end subroutine matrix_read_matrix_market

  !=============================================================================
  subroutine matrix_free_one_shot (f_matrix)
    implicit none
    type(matrix_t), intent(inout) :: f_matrix

    call matrix_free_progressively ( f_matrix, free_only_values )
    call matrix_free_progressively ( f_matrix, free_only_struct )
    call matrix_free_progressively ( f_matrix, free_clean )

  end subroutine matrix_free_one_shot

    !=============================================================================
  subroutine matrix_free_progressively (f_matrix, mode)
    implicit none
    type(matrix_t), intent(inout) :: f_matrix
    integer(ip)     , intent(in)    :: mode 

    if ( mode == free_clean ) then

    else if ( mode == free_only_struct ) then
       ! AFM: In my opinion, the following sentence is 
       ! not convenient at all, at least as a 
       ! general/by-default rule. One might
       ! like to reuse the structure of the graph for
       ! future needs.  Besides, for those cases
       ! where the pointer has been allocated, this would
       ! result in an ugly MEMORY LEAK. I will comment this line 
       ! as a temporal solution, although a cleaner solution
       ! is clearly required.

       ! AFM: the comment above NO longer applies. A sentence
       ! such as allocate(f_matrix%gr) is now FORBIDDEN. The %gr member
       ! will always point to an existing graph, cannot be used to 
       ! allocate a new graph !!! 
       nullify( f_matrix%gr )
    else if ( mode == free_only_values ) then
       call memfree( f_matrix%a,__FILE__,__LINE__)
    end if

  end subroutine matrix_free_progressively

  !=============================================================================
  subroutine matrix_info ( f_mat, me, np )
    implicit none

    ! Parameters 
    type(matrix_t)     , intent(in)    :: f_mat
    integer              , intent(out)   :: me
    integer              , intent(out)   :: np

    me = 0
    np = 1 
  end subroutine matrix_info

  !=============================================================================
  subroutine matrix_zero(mat)
    implicit none
    type(matrix_t), intent(inout)       :: mat

    mat%a = 0.0_rp

  end subroutine matrix_zero

  !=============================================================================
  ! Subroutine to compute C <- A + alpha*B
  subroutine matrix_sum(A,B,C,alpha)
    implicit none
    type(matrix_t)  , intent(in)     :: A,B
    real(rp), optional, intent(in)     :: alpha
    type(matrix_t)  , intent(inout)  :: C

    assert ( A%type == B%type .and. A%type == C%type )
    assert ( A%symm == B%symm .and. A%symm == C%symm )

    if(present(alpha)) then
       C%a = A%a + alpha*B%a
    else
       C%a = A%a + B%a
    end if

  end subroutine matrix_sum

  subroutine matrix_transpose(A, A_t)

    type(matrix_t), intent(in)     :: A         ! Input matrix
    type(matrix_t), intent(out)    :: A_t       ! Output matrix


    ! Locals 
    type(graph_t) :: aux_graph
    integer :: k,i,j

    if ( A%gr%type == csr ) then
       call matrix_alloc ( csr_mat, symm_false,                  &
            A%gr, A_t )
    else if (A%gr%type == csr_symm) then
       call matrix_alloc ( csr_mat, symm_true,                  &
            A%gr, A_t )
    end if

    aux_graph = A%gr

    if (A%gr%type == csr_symm) then 
       A_t%a(:) = A%a(:)
    elseif (A%gr%type == csr ) then
       k = 0    
       do i = 1, A_t%gr%nv
          do j=1, (A_t%gr%ia(i+1) - A_t%gr%ia(i) )
             k = k+1
             A_t%a(k) = A%a(aux_graph%ia( A_t%gr%ja(k) ) )
             aux_graph%ia(A_t%gr%ja(k)) = aux_graph%ia(A_t%gr%ja(k)) + 1
          end do
       end do
    end if

  end subroutine matrix_transpose

  subroutine matvec (a,x,y)
    implicit none
    type(matrix_t) , intent(in)    :: a
    type(vector_t) , intent(in)    :: x
    type(vector_t) , intent(inout) :: y
    real(rp) :: aux

    if (a%symm == symm_false) then
       call matvec_csr(a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
    else if (a%symm == symm_true) then
       call matvec_csr_symm(a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)          
    end if

  end subroutine matvec

  subroutine matmat (a, n, ldX, x, ldY, y)
    implicit none

    ! Parameters
    type(matrix_t) , intent(in)    :: a
    integer(ip)      , intent(in)    :: n
    integer(ip)      , intent(in)    :: ldX
    real(rp)         , intent(in)    :: x(ldX, n)
    integer(ip)      , intent(in)    :: ldY
    real(rp)         , intent(inout) :: y(ldY, n)

    ! Locals 
    integer (ip) :: i

    do i=1,n
       if (a%symm == symm_false) then
          call matvec_csr(a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv2,i),y(1:a%gr%nv,i))
       else if (a%symm == symm_true) then
          call matvec_csr_symm(a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv2,i),y(1:a%gr%nv,i))          
       end if
    end do



  end subroutine matmat

  subroutine matvec_trans (a,x,y)
    implicit none
    type(matrix_t) , intent(in)    :: a
    type(vector_t) , intent(in)    :: x
    type(vector_t) , intent(inout) :: y

    if (a%symm == symm_false) then
       call matvec_csr_trans(a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
    else if (a%symm == symm_true) then
       call matvec_csr_symm_trans(a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)          
    end if

  end subroutine matvec_trans

  subroutine matmat_trans (a, n, ldX, x, ldY, y)
    implicit none

    ! Parameters
    type(matrix_t) , intent(in)    :: a
    integer(ip)      , intent(in)    :: n
    integer(ip)      , intent(in)    :: ldX
    real(rp)         , intent(in)    :: x(ldX, n)
    integer(ip)      , intent(in)    :: ldY
    real(rp)         , intent(inout) :: y(ldY, n)

    ! Locals 
    integer (ip) :: i


    do i=1,n
          if (a%symm == symm_false) then
             call matvec_csr_trans(a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv2,i),y(1:a%gr%nv,i) )
          else if (a%symm == symm_true) then
             call matvec_csr_symm_trans(a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv,i),y(1:a%gr%nv2,i))          
          end if
    end do

  end subroutine matmat_trans

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine matrix_apply(op,x,y) 
    implicit none
    class(matrix_t), intent(in)    :: op
    class(base_operand_t) , intent(in)    :: x
    class(base_operand_t) , intent(inout) :: y 

    call x%GuardTemp()

    select type(x)
    class is (vector_t)
       select type(y)
       class is(vector_t)
          call matvec(op, x, y)
          ! call vector_print(6,y)
       class default
          write(0,'(a)') 'matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'matrix_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()
  end subroutine matrix_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function matrix_apply_fun(op,x) result(y)
    implicit none
    class(matrix_t), intent(in)  :: op
    class(base_operand_t) , intent(in)  :: x
    class(base_operand_t) , allocatable :: y 

    type(vector_t), allocatable :: local_y

    select type(x)
    class is (vector_t)
       allocate(local_y)
       call vector_alloc ( op%gr%nv, local_y)
       call matvec(op, x, local_y)
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'matrix_t%apply_fun: unsupported x class'
       check(1==0)
    end select
  end function matrix_apply_fun

  subroutine matrix_free_tbp(this)
    implicit none
    class(matrix_t), intent(inout) :: this
  end subroutine matrix_free_tbp

end module matrix_names
