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
module serial_scalar_matrix_names
  use types_names
  use memor_names
  use sort_names
  use graph_names
  use serial_scalar_array_names
  
  ! Abstract types
  use vector_names
  use operator_names
  use matrix_names
  use vector_space_names


#ifdef memcheck
  use iso_c_binding
#endif

  implicit none
# include "debug.i90"
  private
  

  
  ! Constants:
  ! Matrix sign
  integer(ip), parameter :: positive_definite     = 0
  integer(ip), parameter :: positive_semidefinite = 1
  integer(ip), parameter :: indefinite            = 2 ! Both positive and negative eigenvalues
  integer(ip), parameter :: unknown               = 3 ! No info
  
  type, extends(matrix_t) :: serial_scalar_matrix_t
     logical                    :: is_symmetric    
     integer(ip)                :: sign            
     real(rp)     , allocatable :: a(:)            
     type(graph_t), pointer     :: graph => NULL() 
   contains
     procedure, private :: serial_scalar_matrix_create_square
     procedure, private :: serial_scalar_matrix_create_rectangular
     generic  :: create => serial_scalar_matrix_create_square, &
          & serial_scalar_matrix_create_rectangular

     procedure  :: allocate                        => serial_scalar_matrix_allocate
     procedure  :: print                           => serial_scalar_matrix_print
     procedure  :: print_matrix_market             => serial_scalar_matrix_print_matrix_market
     procedure  :: read_matrix_market              => serial_scalar_matrix_read_matrix_market
     procedure  :: apply_to_dense_matrix           => serial_scalar_matrix_apply_to_dense_matrix
     procedure  :: apply_transpose_to_dense_matrix => serial_scalar_matrix_apply_transpose_to_dense_matrix 
     procedure  :: transpose                       => serial_scalar_matrix_transpose
     procedure  :: init                            => serial_scalar_matrix_init
     procedure  :: apply                           => serial_scalar_matrix_apply
     procedure  :: free_in_stages                  => serial_scalar_matrix_free_in_stages
     procedure  :: default_initialization          => serial_scalar_matrix_default_init
  end type serial_scalar_matrix_t

  ! Constants
  public :: positive_definite, positive_semidefinite, indefinite, unknown

  ! Types
  public :: serial_scalar_matrix_t

!***********************************************************************
! Allocatable arrays of type(serial_scalar_matrix_t)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(serial_scalar_matrix_t)
# define var_size 80
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine serial_scalar_matrix_create_square(this,symmetric_storage,is_symmetric,sign)
    implicit none
    class(serial_scalar_matrix_t), intent(out) :: this
    logical                     , intent(in)  :: symmetric_storage
    logical                     , intent(in)  :: is_symmetric
    integer(ip)                 , intent(in)  :: sign

    integer(ip) :: istat

    allocate( this%graph, stat=istat )
    check ( istat == 0 )
    call this%graph%create ( symmetric_storage )	

    assert ( sign == positive_definite .or. sign == positive_semidefinite  .or. sign == indefinite .or. sign == unknown )
    this%is_symmetric = is_symmetric
    this%sign = sign
  end subroutine serial_scalar_matrix_create_square

  !=============================================================================
  subroutine serial_scalar_matrix_create_rectangular(this)
    implicit none
    class(serial_scalar_matrix_t), intent(out) :: this
    integer(ip) :: istat
    allocate( this%graph, stat=istat )
    check ( istat == 0 )
    call this%graph%create ( symmetric_storage=.false. )	
    this%is_symmetric = .false.
    this%sign = unknown
  end subroutine serial_scalar_matrix_create_rectangular

  subroutine serial_scalar_matrix_allocate(this)
    implicit none
    class(serial_scalar_matrix_t), intent(inout) :: this
    type(serial_scalar_array_t) :: range_vector
    type(serial_scalar_array_t) :: domain_vector
    type(vector_space_t), pointer :: range_vector_space
    type(vector_space_t), pointer :: domain_vector_space

    call memalloc(this%graph%ia(this%graph%nv+1)-1,this%a,__FILE__,__LINE__)
    this%a = 0.0_rp
    
    call range_vector%create(this%graph%nv)
    call domain_vector%create(this%graph%nv2)

    range_vector_space => this%get_range_vector_space()
    call range_vector_space%create(range_vector)
    
    domain_vector_space => this%get_domain_vector_space()
    call domain_vector_space%create(domain_vector)
    
    call range_vector%free()
    call domain_vector%free()
  end subroutine serial_scalar_matrix_allocate

  !=============================================================================
  subroutine serial_scalar_matrix_free_in_stages (this, action)
    implicit none
    class(serial_scalar_matrix_t), intent(inout) :: this
    integer(ip)                  , intent(in)    :: action
    integer(ip) :: istat
    if ( action == free_clean ) then
       deallocate( this%graph, stat=istat )
       check ( istat == 0 )
    else if ( action == free_struct ) then
       call this%graph%free() 
    else if ( action == free_values ) then
       call memfree( this%a,__FILE__,__LINE__)
       call this%free_vector_spaces()
    end if
  end subroutine serial_scalar_matrix_free_in_stages

  !=============================================================================
  subroutine serial_scalar_matrix_default_init (this)
    implicit none
    class(serial_scalar_matrix_t), intent(inout) :: this
    this%is_symmetric = .false.
    this%sign = unknown
    nullify(this%graph)
    call this%NullifyTemporary()
  end subroutine serial_scalar_matrix_default_init

  !=============================================================================
  subroutine serial_scalar_matrix_print(this,lunou)
    implicit none
    class(serial_scalar_matrix_t), intent(in)    :: this
    integer(ip)     , intent(in)    :: lunou
    integer(ip)                     :: i
    assert ( associated(this%graph) )
    call this%graph%print (lunou)
    write (lunou, '(a)')     '*** begin matrix data structure ***'
    do i=1,this%graph%nv
       write(lunou,'(10(e25.16,1x))') this%a(this%graph%ia(i):this%graph%ia(i+1)-1)
    end do
  end subroutine serial_scalar_matrix_print

  subroutine serial_scalar_matrix_print_matrix_market (this, lunou, ng, l2g)
    implicit none
    class(serial_scalar_matrix_t), intent(in)           :: this
    integer(ip)     , intent(in)           :: lunou
    integer(ip)     , intent(in), optional :: ng
    integer(ip)     , intent(in), optional :: l2g (*)

    integer(ip) :: i, j
    integer(ip) :: nv1_, nv2_, nl_


    if ( present(ng) ) then 
       nv1_ = ng
       nv2_ = ng
    else
       nv1_ = this%graph%nv
       nv2_ = this%graph%nv2
    end if

    write (lunou,'(a)') '%%MatrixMarket matrix coordinate real general'
    if (.not. this%graph%symmetric_storage) then
       write (lunou,*) nv1_,nv2_,this%graph%ia(this%graph%nv+1)-1
       do i=1,this%graph%nv
          do j=this%graph%ia(i),this%graph%ia(i+1)-1
             if (present(l2g)) then
                write(lunou,'(i12, i12, e32.25)') l2g(i), l2g(this%graph%ja(j)), this%a(j)
             else
                write(lunou,'(i12, i12, e32.25)') i, this%graph%ja(j), this%a(j)
             end if
          end do
       end do
    else 
       write (lunou,*) nv1_,nv2_,& 
            2*(this%graph%ia(this%graph%nv+1)-1) - this%graph%nv

       do i=1,this%graph%nv
          do j=this%graph%ia(i),this%graph%ia(i+1)-1
             if (present(l2g)) then
                write(lunou,'(i12, i12, e32.25)') l2g(i), l2g(this%graph%ja(j)), this%a(j)
             else
                write(lunou,'(i12, i12, e32.25)') i, this%graph%ja(j), this%a(j)
             end if
             if (i /= this%graph%ja(j)) then
                if (present(l2g)) then
                   write(lunou,'(i12, i12, e32.25)') l2g(this%graph%ja(j)), l2g(i), this%a(j)
                else
                   write(lunou,'(i12, i12, e32.25)') this%graph%ja(j), i, this%a(j)
                end if
             end if
          end do
       end do
    end if

  end subroutine serial_scalar_matrix_print_matrix_market

  subroutine serial_scalar_matrix_read_matrix_market (this, lunou, symmetric_storage, is_symmetric, sign)
    implicit none
    class(serial_scalar_matrix_t), intent(out)     :: this
    integer(ip)                  , intent(in)      :: lunou
    logical                      , intent(in)      :: symmetric_storage
    logical                      , intent(in)      :: is_symmetric
    integer(ip)                  , intent(in)      :: sign

    integer(ip) :: i, j, iw1(2),iw2(2)
    integer(ip) :: nv, nw, nz

    integer(ip), allocatable :: ija_work(:,:), ija_index(:)
    real(rp)   , allocatable :: a_work(:)

    call this%create(symmetric_storage, is_symmetric, sign)

    read (lunou,*) 
    read (lunou,*) nv,nw,nz

    if(nv/=nw) then 
       write (0,*) 'Error: only square matrices are allowed in fempar when reading from matrix market'
       stop
    end if

    call memalloc ( 2, nz, ija_work, __FILE__,__LINE__ )
    call memalloc (    nz, ija_index, __FILE__,__LINE__ )
    call memalloc (    nz, a_work,__FILE__,__LINE__)

    if(.not. this%graph%symmetric_storage) then
       do i=1,nz
          read(lunou,'(i12, i12, e32.25)',end=10, err=10) ija_work(1,i), ija_work(2,i), a_work(i)
          ija_index(i)=i
       end do
    else
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
    this%graph%nv   = nv
    this%graph%nv2  = nv
    call memalloc ( nv+1 , this%graph%ia, __FILE__,__LINE__ )
    call memalloc ( nz , this%graph%ja, __FILE__,__LINE__ )

    ! Count columns on each row
    this%graph%ia = 0
    do i=1,nz
       this%graph%ia(ija_work(1,i)+1) = this%graph%ia(ija_work(1,i)+1) + 1
    end do
    ! Compress
    this%graph%ia(1)=1
    do i=1,nv
       this%graph%ia(i+1) = this%graph%ia(i+1) +this%graph%ia(i)
    end do
    write(*,*) this%graph%ia(nv+1)
    ! Copy ja
    this%graph%ja = ija_work(2,:)
    call memfree(ija_work,__FILE__,__LINE__)

    ! Reorder a
    call this%allocate()
    do i=1,nz
       this%a(i)=a_work(ija_index(i))
    end do
    call memfree ( ija_index,__FILE__,__LINE__)
    call memfree (    a_work,__FILE__,__LINE__)

    return

10  write (0,*) 'Error reading matrix eof or err'
    check(.false.)

  end subroutine serial_scalar_matrix_read_matrix_market

  !=============================================================================
  subroutine serial_scalar_matrix_init(this, alpha)
    implicit none
    class(serial_scalar_matrix_t), intent(inout) :: this
    real(rp)                     , intent(in)    :: alpha
    this%a = alpha
  end subroutine serial_scalar_matrix_init

  subroutine serial_scalar_matrix_transpose(this, A_t)
    class(serial_scalar_matrix_t), intent(in)     :: this    ! Input matrix
    type(serial_scalar_matrix_t), intent(inout)  :: A_t  ! Output matrix

    ! Locals 
    type(graph_t) :: aux_graph
    integer :: k,i,j

    assert ( this%graph%symmetric_storage .eqv. A_t%graph%symmetric_storage )
    call this%graph%copy(A_t%graph)
    call A_t%allocate()

    if (this%graph%symmetric_storage) then 
       A_t%a(:) = this%a(:)
    else
       call this%graph%copy(aux_graph )
       k = 0    
       do i = 1, A_t%graph%nv
          do j=1, (A_t%graph%ia(i+1) - A_t%graph%ia(i) )
             k = k+1
             A_t%a(k) = this%a(aux_graph%ia( A_t%graph%ja(k) ) )
             aux_graph%ia(A_t%graph%ja(k)) = aux_graph%ia(A_t%graph%ja(k)) + 1
          end do
       end do
       call aux_graph%free()
    end if
  end subroutine serial_scalar_matrix_transpose

  subroutine serial_scalar_matrix_matvec (a,x,y)
    implicit none
    type(serial_scalar_matrix_t) , intent(in)    :: a
    type(serial_scalar_array_t) , intent(in)    :: x
    type(serial_scalar_array_t) , intent(inout) :: y
    real(rp) :: aux

    if (.not. a%graph%symmetric_storage) then
       call matvec(a%graph%nv,a%graph%nv2,a%graph%ia,a%graph%ja,a%a,x%b,y%b)
    else 
       call matvec_symmetric_storage(a%graph%nv,a%graph%nv,a%graph%ia,a%graph%ja,a%a,x%b,y%b)          
    end if

  end subroutine serial_scalar_matrix_matvec

  subroutine serial_scalar_matrix_apply_to_dense_matrix (this, n, ldX, x, ldY, y)
    implicit none
    class(serial_scalar_matrix_t) , intent(in)    :: this
    integer(ip)                   , intent(in)    :: n
    integer(ip)                   , intent(in)    :: ldX
    real(rp)                      , intent(in)    :: x(ldX, n)
    integer(ip)                   , intent(in)    :: ldY
    real(rp)                      , intent(inout) :: y(ldY, n)

    ! Locals 
    integer (ip) :: i

    do i=1,n
       if (.not. this%graph%symmetric_storage) then
          call matvec(this%graph%nv,this%graph%nv2,this%graph%ia,this%graph%ja,this%a,x(1:this%graph%nv2,i),y(1:this%graph%nv,i))
       else 
          call matvec_symmetric_storage(this%graph%nv,this%graph%nv,this%graph%ia,this%graph%ja,this%a,x(1:this%graph%nv2,i),y(1:this%graph%nv,i))          
       end if
    end do
  end subroutine serial_scalar_matrix_apply_to_dense_matrix

  subroutine serial_scalar_matrix_apply_transpose_to_dense_matrix (this, n, ldX, x, ldY, y)
    implicit none

    ! Parameters
    class(serial_scalar_matrix_t) , intent(in)    :: this
    integer(ip)                   , intent(in)    :: n
    integer(ip)                   , intent(in)    :: ldX
    real(rp)                      , intent(in)    :: x(ldX, n)
    integer(ip)                   , intent(in)    :: ldY
    real(rp)                      , intent(inout) :: y(ldY, n)

    ! Locals 
    integer (ip) :: i

    do i=1,n
       if (.not. this%graph%symmetric_storage) then
          call matvec_trans(this%graph%nv,this%graph%nv2,this%graph%ia,this%graph%ja,this%a,x(1:this%graph%nv2,i),y(1:this%graph%nv,i) )
       else 
          call matvec_symmetric_storage_trans(this%graph%nv,this%graph%nv,this%graph%ia,this%graph%ja,this%a,x(1:this%graph%nv,i),y(1:this%graph%nv2,i))          
       end if
    end do

  end subroutine serial_scalar_matrix_apply_transpose_to_dense_matrix

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine serial_scalar_matrix_apply(op,x,y) 
    implicit none
    class(serial_scalar_matrix_t), intent(in)    :: op
    class(vector_t) , intent(in)    :: x
    class(vector_t) , intent(inout) :: y 

    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    
    call x%GuardTemp()
    select type(x)
       class is (serial_scalar_array_t)
       select type(y)
          class is(serial_scalar_array_t)
          call serial_scalar_matrix_matvec(op, x, y)
       end select
    end select
    call x%CleanTemp()
  end subroutine serial_scalar_matrix_apply

  ! Debugged
  subroutine matvec (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
    real(rp)   , intent(out) :: y(nv)
    integer(ip)              :: iv,iz,jv

    y = 0.0_rp
    do iv = 1, nv
       do iz = ia(iv), ia(iv+1)-1
          jv   = ja(iz)
          y(iv) = y(iv) + x(jv)*a(iz)
       end do ! iz
    end do ! iv

  end subroutine matvec

  ! Debugged
  subroutine matvec_symmetric_storage (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
    real(rp)   , intent(out) :: y(nv)
    integer(ip)              :: iv,iz,jv

    assert(nv==nv2)

    y = 0.0_rp
    do iv = 1, nv
       y(iv) = y(iv) + x(ja(ia(iv)))*a(ia(iv))
       do iz = ia(iv)+1, ia(iv+1)-1
          jv = ja(iz)
          y(iv) = y(iv) + x(jv)*a(iz)
          y(jv) = y(jv) + x(iv)*a(iz)
       end do ! iz
    end do ! iv

  end subroutine matvec_symmetric_storage

  ! Debugged
  subroutine matvec_trans (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
    real(rp)   , intent(out) :: y(nv2)
    integer(ip)              :: iv,iz,jv

    y = 0.0_rp
    do iv = 1, nv
       do iz = ia(iv), ia(iv+1)-1
          jv = ja(iz)
          y(jv) = y(jv) + x(iv)*a(iz)
       end do ! iz
    end do ! iv

  end subroutine matvec_trans

  subroutine matvec_symmetric_storage_trans (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
    real(rp)   , intent(out) :: y(nv2)
    integer(ip)              :: iv,iz,jv,id,jd
    integer(ip)              :: of, ivc, izc, jvc

    write (0,*) 'Error: the body of matvec_symmetric_storage_trans in matvec.f90 still to be written'
    write (0,*) 'Error: volunteers are welcome !!!'
    stop 
  end subroutine matvec_symmetric_storage_trans

end module serial_scalar_matrix_names
