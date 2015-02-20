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
module fem_matrix_class
  use types
  use memor
  use sort_class
  use fem_graph_class
  use elmat_class
  use matrix_dof
  use dof_handler_class
  use fem_space_class
  use fem_blocks_class
  use adaptivity
#ifdef memcheck
  use iso_c_binding
#endif
  ! CRUCIAL NOTE: matrix requires BOUNDARY conditions
  ! temporarily for bnn/bddc dd algorithms.
  ! A more clean/definitive solution is required
  ! here to treat these situations.
  use fem_conditions_class

  ! CRUCIAL NOTE: matrix requires BOUNDARY conditions
  ! temporarily for bnn/bddc dd algorithms.
  ! A more clean/definitive solution is required
  ! here to treat these situations.
  use fem_mesh_class

  implicit none
# include "debug.i90"

  private

  ! Constants:
  ! a) fem_matrix types
  integer(ip), parameter :: css_mat=10
  integer(ip), parameter :: csr_mat=20
  integer(ip), parameter :: csc_mat=30
  ! b) fem_matrix symmetry
  integer(ip), parameter :: symm_true =0
  integer(ip), parameter :: symm_false=1
  ! c) fem_matrix sign
  integer(ip), parameter :: positive_definite     = 0
  integer(ip), parameter :: positive_semidefinite = 1
  integer(ip), parameter :: indefinite            = 2 ! Both positive and negative eigenvalues
  integer(ip), parameter :: unknown               = 3 ! No info

  ! Matrix
  type fem_matrix
     integer(ip)                :: &
          storage=undef_sto,         &         ! Storage layout (blk: block; scal: scalar)
          symm=symm_false,           &         ! Flag for symmetry
          sign=positive_definite,    &         ! Flag for positiveness
          type=css_mat,              &         ! fem_matrix type (csr, csc, css, epetra, petsc )
          nd1=0,                     &         ! Number of degrees of freedom, ndof1
          nd2=0                                ! Number of degrees of freedom, ndof2


     ! We need to decide how blocks will be stored (transposed or not)
     ! Currently in felap u is transposed but l is not. Seems better to
     ! transpose both. What about trilinos?
     real(rp), allocatable      :: &
          d(:,:,:),                  &         ! Diagonal components (nd1,nd2,neq) 
          l(:,:,:),                  &         ! Lower in css (nd1,nd2,nzs)        
          u(:,:,:),                  &         ! Upper in css (nd1,nd2,nzs)        
          a(:,:,:)                             ! Lower/Upper part components ordered as:
     ! Rows (nd1,nd2,nzt)                
     ! Cols (nd2,nd1,nzt)                

     !type(epetra_csrmat) :: emp

     type(fem_graph),pointer :: &
          gr => NULL()                      ! Associated fem_graph

     ! CRUCIAL NOTE: matrix requires BOUNDARY conditions
     ! temporarily for bnn/bddc dd algorithms.
     ! A more clean/definitive solution is required
     ! here to treat these situations.
     type(fem_conditions),pointer :: &
          bcs => NULL()                      ! Associated fem_conditions

     ! CRUCIAL NOTE: matrix requires mesh point coordinates
     ! temporarily for bnn/bddc dd algorithms.
     ! A more clean/definitive solution is required
     ! here to treat these situations.
     type(fem_mesh),pointer :: &
          msh => NULL()                      ! Associated fem_mesh

  end type fem_matrix

  interface fem_matrix_assembly
     module procedure fem_matrix_assembly_standard
     module procedure fem_matrix_assembly_w_dof_handler
     module procedure fem_matrix_assembly_w_dof_coupl
     module procedure fem_matrix_assembly_nosq_w_dof_handler
  end interface fem_matrix_assembly

  ! Memalloc interfaces
  interface memalloc
     module procedure memalloc_fem_matrix1, memalloc_fem_matrix2, memalloc_fem_matrix3, &
          &           memalloc_fem_matrix4
  end interface memalloc
  interface memrealloc
     module procedure memrealloc_fem_matrix1, memrealloc_fem_matrix2, memrealloc_fem_matrix3, &
          &           memrealloc_fem_matrix4
  end interface memrealloc
  interface memfree
     module procedure memfree_fem_matrix1, memfree_fem_matrix2, memfree_fem_matrix3, &
          &           memfree_fem_matrix4
  end interface memfree
  interface memmovealloc
     module procedure memmovealloc_fem_matrix1, memmovealloc_fem_matrix2, memmovealloc_fem_matrix3, &
          &           memmovealloc_fem_matrix4
  end interface memmovealloc

  interface fem_matrix_free
     module procedure fem_matrix_free_one_shot, fem_matrix_free_progressively
  end interface fem_matrix_free

  ! Constants
  public :: css_mat, csr_mat, csc_mat, symm_true, symm_false
  public :: positive_definite, positive_semidefinite, indefinite, unknown

  ! Types
  public :: fem_matrix

  ! Functions
  public :: fem_matrix_create, fem_matrix_graph, fem_matrix_fill_val, &
       &    fem_matrix_alloc, fem_matrix_copy, fem_matrix_print, fem_matrix_read,    &
       &    fem_matrix_set_bcs, fem_matrix_set_msh, fem_matrix_print_matrix_market,  & 
       &    fem_matrix_assembly, fem_matrix_free, fem_matrix_info, fem_matrix_zero,  &
       &    fem_matrix_sum, fem_matrix_transpose, &
       &    fem_matrix_compose_name_matrix_market, fem_matrix_read_matrix_market,    &
       &    fem_face_matrix_assembly_w_dof_handler, fem_matrix_assembly_w_dof_coupl, &
       &    memalloc, memrealloc, memfree, memmovealloc

contains

  !***********************************************************************
  !***********************************************************************
  ! Specialization to allocate type(fem_matrix)
  !***********************************************************************
  !***********************************************************************
# define var_type type(fem_matrix)
# define var_size 80
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_fem_matrix1    
# define generic_memalloc_2      memalloc_fem_matrix2    
# define generic_memalloc_3      memalloc_fem_matrix3    
# define generic_memalloc_4      memalloc_fem_matrix4    
# define generic_memrealloc_1    memrealloc_fem_matrix1  
# define generic_memrealloc_2    memrealloc_fem_matrix2  
# define generic_memrealloc_3    memrealloc_fem_matrix3  
# define generic_memrealloc_4    memrealloc_fem_matrix4  
# define generic_memfree_1       memfree_fem_matrix1     
# define generic_memfree_2       memfree_fem_matrix2     
# define generic_memfree_3       memfree_fem_matrix3     
# define generic_memfree_4       memfree_fem_matrix4     
# define generic_memmovealloc_1  memmovealloc_fem_matrix1
# define generic_memmovealloc_2  memmovealloc_fem_matrix2
# define generic_memmovealloc_3  memmovealloc_fem_matrix3
# define generic_memmovealloc_4  memmovealloc_fem_matrix4
# include "memor.i90"

  !=============================================================================
  subroutine fem_matrix_create(storage,type,symm,nd1,nd2,mat,def)
    implicit none
    integer(ip)     , intent(in)           :: storage, type, symm, nd1, nd2
    type(fem_matrix), intent(out)          :: mat
    integer(ip)     , optional, intent(in) :: def

    assert ( storage == blk    .or. storage == scal )
    assert ( type == css_mat   .or. type == csr_mat .or. type == csc_mat )
    assert ( symm == symm_true .or. symm == symm_false )
    assert ( nd1 >= 1 .and. nd2 >= 1 ) 

    mat%storage =  storage ! Storage layout (blk: block; scal: scalar)
    mat%symm    =  symm    ! Flag for symmetry
    mat%type    =  type    ! fem_matrix type (csr_mat, csc_mat, css_mat)
    mat%nd1     =  nd1     ! Number of degrees of freedom, ndof1
    mat%nd2     =  nd2     ! Number of degrees of freedom, ndof2

    mat%sign = unknown
    if(present(def)) then
       assert ( def == positive_definite .or. def == positive_semidefinite  .or. def == indefinite .or. def == unknown )
       mat%sign = def
    end if

  end subroutine fem_matrix_create

  subroutine fem_matrix_graph(gr,mat)
    implicit none
    type(fem_graph) , target, intent(in) :: gr
    type(fem_matrix), intent(inout)      :: mat
    mat%gr => gr
  end subroutine fem_matrix_graph

  subroutine fem_matrix_fill_val(mat)
    implicit none
    type(fem_matrix), intent(inout)      :: mat
    integer(ip)                          :: neq, nzs, nzt

    neq  =  mat%gr%nv              ! Number of rows and columns (equations)

    if ( mat%storage == blk ) then ! gr must describe block sparsity pattern of mat
       if(mat%type==css_mat) then
          assert(mat%nd1==mat%nd2)
          nzs=mat%gr%nzs
          call memalloc(mat%nd1,mat%nd2,neq,mat%d,__FILE__,__LINE__)
          mat%d = 0.0_rp
          call memalloc(mat%nd1,mat%nd2,nzs,mat%l,__FILE__,__LINE__)
          mat%l = 0.0_rp
          if(mat%symm==symm_true) then
             ! only an empty address to have same interfaces
             call memalloc(1,1,1,mat%u,__FILE__,__LINE__)
          else if(mat%symm==symm_false) then
             call memalloc(mat%nd2,mat%nd1,nzs,mat%u,__FILE__,__LINE__)
             mat%u = 0.0_rp
          end if
       else if(mat%type==csr_mat) then
          nzt = mat%gr%ia(mat%gr%nv+1)-1
          call memalloc(mat%nd2,mat%nd1,nzt,mat%a,__FILE__,__LINE__)
          mat%a = 0.0_rp
       else if(mat%type==csc_mat) then
          nzt = mat%gr%ia(mat%gr%nv+1)-1
          call memalloc(mat%nd1,mat%nd2,nzt,mat%a,__FILE__,__LINE__)
          mat%a = 0.0_rp
       end if
    else if (mat%storage == scal) then ! gr must describe scalar sparsity pattern of mat
       if(mat%type==css_mat) then
          assert(mat%nd1==mat%nd2)
          nzs=mat%gr%nzs
          call memalloc(1,1,neq,mat%d,__FILE__,__LINE__)
          mat%d = 0.0_rp
          call memalloc(1,1,nzs,mat%l,__FILE__,__LINE__)
          mat%l = 0.0_rp
          if(mat%symm==symm_true) then
             ! only an empty address to have same interfaces
             call memalloc(1,1,1,mat%u,__FILE__,__LINE__)
          else if(mat%symm==symm_false) then
             call memalloc(1,1,nzs,mat%u,__FILE__,__LINE__)
             mat%u = 0.0_rp
          end if
       else if(mat%type==csr_mat) then
          nzt = mat%gr%ia(mat%gr%nv+1)-1
          call memalloc(1,1,nzt,mat%a,__FILE__,__LINE__)
          mat%a = 0.0_rp
       else if(mat%type==csc_mat) then
          nzt = mat%gr%ia(mat%gr%nv+1)-1
          call memalloc(1,1,nzt,mat%a,__FILE__,__LINE__)
          mat%a = 0.0_rp
       end if
    end if
  end subroutine fem_matrix_fill_val

  subroutine fem_matrix_alloc(storage,type,symm,nd1,nd2,gr,mat,def)
    implicit none
    integer(ip)     , intent(in)           :: storage, type, symm, nd1, nd2
    type(fem_graph) , target, intent(in)   :: gr
    type(fem_matrix), intent(out)          :: mat
    integer(ip)     , optional, intent(in) :: def

    call fem_matrix_create(storage,type,symm,nd1,nd2,mat,def)
    call fem_matrix_graph(gr,mat)
    call fem_matrix_fill_val(mat)

    ! integer(ip)                           :: neq, nzs, nzt

    ! assert ( storage == blk    .or. storage == scal )
    ! assert ( type == css_mat   .or. type == csr_mat .or. type == csc_mat )
    ! assert ( symm == symm_true .or. symm == symm_false )
    ! assert ( nd1 >= 1 .and. nd2 >= 1 ) 

    ! mat%storage =  storage ! Storage layout (blk: block; scal: scalar)
    ! mat%symm    =  symm    ! Flag for symmetry
    ! mat%type    =  type    ! fem_matrix type (csr_mat, csc_mat, css_mat)
    ! mat%nd1     =  nd1     ! Number of degrees of freedom, ndof1
    ! mat%nd2     =  nd2     ! Number of degrees of freedom, ndof2
    ! neq         =  gr%nv   ! Number of rows and columns (equations)

    ! mat%gr      => gr

    ! mat%sign = unknown
    ! if(present(def)) then
    !    assert ( def == positive_definite .or. def == positive_semidefinite  .or. def == indefinite .or. def == unknown )
    !    mat%sign = def
    ! end if

    ! if ( storage == blk ) then ! gr must describe block sparsity pattern of mat
    !    if(type==css_mat) then
    !       assert(mat%nd1==mat%nd2)
    !       nzs=gr%nzs
    !       call memalloc(mat%nd1,mat%nd2,neq,mat%d,__FILE__,__LINE__)
    !       mat%d = 0.0_rp
    !       call memalloc(mat%nd1,mat%nd2,nzs,mat%l,__FILE__,__LINE__)
    !       mat%l = 0.0_rp
    !       if(mat%symm==symm_true) then
    !          ! only an empty address to have same interfaces
    !          call memalloc(1,1,1,mat%u,__FILE__,__LINE__)
    !       else if(mat%symm==symm_false) then
    !          call memalloc(mat%nd2,mat%nd1,nzs,mat%u,__FILE__,__LINE__)
    !          mat%u = 0.0_rp
    !       end if
    !    else if(type==csr_mat) then
    !       nzt = mat%gr%ia(mat%gr%nv+1)-1
    !       call memalloc(mat%nd2,mat%nd1,nzt,mat%a,__FILE__,__LINE__)
    !       mat%a = 0.0_rp
    !    else if(type==csc_mat) then
    !       nzt = mat%gr%ia(mat%gr%nv+1)-1
    !       call memalloc(mat%nd1,mat%nd2,nzt,mat%a,__FILE__,__LINE__)
    !       mat%a = 0.0_rp
    !    end if
    ! else if (storage == scal) then ! gr must describe scalar sparsity pattern of mat
    !    if(type==css_mat) then
    !       assert(mat%nd1==mat%nd2)
    !       nzs=gr%nzs
    !       call memalloc(1,1,neq,mat%d,__FILE__,__LINE__)
    !       mat%d = 0.0_rp
    !       call memalloc(1,1,nzs,mat%l,__FILE__,__LINE__)
    !       mat%l = 0.0_rp
    !       if(mat%symm==symm_true) then
    !          ! only an empty address to have same interfaces
    !          call memalloc(1,1,1,mat%u,__FILE__,__LINE__)
    !       else if(mat%symm==symm_false) then
    !          call memalloc(1,1,nzs,mat%u,__FILE__,__LINE__)
    !          mat%u = 0.0_rp
    !       end if
    !    else if(type==csr_mat) then
    !       nzt = mat%gr%ia(mat%gr%nv+1)-1
    !       call memalloc(1,1,nzt,mat%a,__FILE__,__LINE__)
    !       mat%a = 0.0_rp
    !    else if(type==csc_mat) then
    !       nzt = mat%gr%ia(mat%gr%nv+1)-1
    !       call memalloc(1,1,nzt,mat%a,__FILE__,__LINE__)
    !       mat%a = 0.0_rp
    !    end if
    ! end if
  end subroutine fem_matrix_alloc

  subroutine fem_matrix_copy (imatrix, omatrix)
    implicit none
    ! Parameters 
    type(fem_matrix), intent(in)    :: imatrix
    type(fem_matrix), intent(inout) :: omatrix

    ! *** IMPORTANT NOTE: This routine assumes that omatrix 
    ! already has an associated graph in omatrix%gr

    ! This routine is only provided for matrices stored in CSR
    ! format. If you need to copy a matrix in a different format,
    ! please extend this routine
    assert ( imatrix%type == csr_mat )

    call fem_matrix_alloc ( imatrix%storage, imatrix%type, imatrix%symm, imatrix%nd1, imatrix%nd2, omatrix%gr, omatrix, imatrix%sign)
    omatrix%a = imatrix%a

  end subroutine fem_matrix_copy


  !=============================================================================
  subroutine fem_matrix_set_bcs (bcs, mat)
    implicit none
    type(fem_conditions), intent(in), target :: bcs
    type(fem_matrix), intent(inout)          :: mat

    ! CRUCIAL NOTE: matrix requires BOUNDARY conditions
    ! temporarily for bnn/bddc dd algorithms.
    ! A more clean/definitive solution is required
    ! here to treat these situations.
    mat%bcs => bcs 
  end subroutine fem_matrix_set_bcs

  !=============================================================================
  subroutine fem_matrix_set_msh (msh, mat)
    implicit none
    type(fem_mesh), intent(in), target :: msh
    type(fem_matrix), intent(inout)    :: mat

    ! CRUCIAL NOTE: matrix requires BOUNDARY conditions
    ! temporarily for bnn/bddc dd algorithms.
    ! A more clean/definitive solution is required
    ! here to treat these situations.
    mat%msh => msh 
  end subroutine fem_matrix_set_msh

  !=============================================================================
  subroutine fem_matrix_print(lunou, f_matrix)
    implicit none
    type(fem_matrix), intent(in)    :: f_matrix
    integer(ip)     , intent(in)    :: lunou
    integer(ip)                     :: i

    assert ( associated(f_matrix%gr) )
    call fem_graph_print (lunou, f_matrix%gr)

    write (lunou, '(a)')     '*** begin fem_matrix data structure ***'
    if ( f_matrix%type == csr_mat ) then
       ! write (lunou, '(4E20.13)')  f_matrix%a(:, :, :)
       do i=1,f_matrix%gr%nv
          !write(lunou,'(10(e14.7,1x))') f_matrix%a(:,:,f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
          write(lunou,'(10(e25.16,1x))') f_matrix%a(:,:,f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
       end do
    end if

  end subroutine fem_matrix_print

  !=============================================================================
  subroutine fem_matrix_read(lunin, f_matrix)
    implicit none
    ! Parameters
    integer(ip)     , intent(in)    :: lunin
    type(fem_matrix), intent(inout) :: f_matrix
    
    ! Locals
    integer(ip)                     :: i,nzt

    ! AFM: the following sentence is now NO longer permitted.
    !      f_matrix should be passed to fem_matrix_read in
    !      such a way that the graph is already associated 
    ! allocate(f_matrix%gr)
    ! call fem_graph_read (lunin, f_matrix%gr)
    assert ( associated(f_matrix%gr) )

    if ( f_matrix%type == csr_mat ) then
       nzt = f_matrix%gr%ia(f_matrix%gr%nv+1)-1
       call memalloc(1,1,nzt,f_matrix%a,__FILE__,__LINE__)

       do i=1,f_matrix%gr%nv/10
          read(lunin,'(10e15.7)') f_matrix%a(1,1,f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
       end do
       if(mod(f_matrix%gr%nv,10)>0) &
            &read(lunin,'(10e15.7)') f_matrix%a(1,1,f_matrix%gr%ia(i):f_matrix%gr%ia(i+1)-1)
    end if

  end subroutine fem_matrix_read

  !=============================================================================
  subroutine fem_matrix_compose_name_matrix_market ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.mtx'
  end subroutine fem_matrix_compose_name_matrix_market

  subroutine fem_matrix_print_matrix_market (lunou, f_matrix, ng, l2g)
    implicit none
    type(fem_matrix), intent(in)           :: f_matrix
    integer(ip)     , intent(in)           :: lunou
    integer(ip)     , intent(in), optional :: ng
    integer(ip)     , intent(in), optional :: l2g (*)

    integer(ip) :: i, j
    integer(ip) :: nv1_, nv2_, nl_

    if (f_matrix%storage == blk) then ! gr must describe block sparsity pattern of f_matrix
       if(f_matrix%type==css_mat) then
          write (0,*) 'Error: the body of fem_matrix_print_matrix_market in matrix.f90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       else if(f_matrix%type==csr_mat) then
          write (0,*) 'Error: the body of fem_matrix_print_matrix_market in matrix.f90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       else if(f_matrix%type==csc_mat) then
          write (0,*) 'Error: the body of fem_matrix_print_matrix_market in matrix.f90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       end if
    else if (f_matrix%storage == scal) then ! gr must describe scalar sparsity pattern of f_matrix
       if( f_matrix%type==css_mat ) then
          write (0,*) 'Error: the body of fem_matrix_print_matrix_market in matrix.f90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       else if(f_matrix%type==csr_mat) then

          ! assert (f_matrix%nd1==1.and.f_matrix%nd2==1)
          assert ( f_matrix%storage == scal )

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
                      write(lunou,'(i10, i10, e32.25)') l2g(i), l2g(f_matrix%gr%ja(j)), f_matrix%a(1,1,j)
                   else
                      write(lunou,'(i10, i10, e32.25)') i, f_matrix%gr%ja(j), f_matrix%a(1,1,j)
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
                      write(lunou,'(i10, i10, e32.25)') l2g(i), l2g(f_matrix%gr%ja(j)), f_matrix%a(1,1,j)
                   else
                      write(lunou,'(i10, i10, e32.25)') i, f_matrix%gr%ja(j), f_matrix%a(1,1,j)
                   end if
                   if (i /= f_matrix%gr%ja(j)) then
                      if (present(l2g)) then
                         write(lunou,'(i10, i10, e32.25)') l2g(f_matrix%gr%ja(j)), l2g(i), f_matrix%a(1,1,j)
                      else
                         write(lunou,'(i10, i10, e32.25)') f_matrix%gr%ja(j), i, f_matrix%a(1,1,j)
                      end if
                   end if
                end do
             end do
          end if
       else if(f_matrix%type==csc_mat) then
          write (0,*) 'Error: the body of fem_matrix_print_matrix_market in matrix.f90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       end if
    end if

  end subroutine fem_matrix_print_matrix_market

  subroutine fem_matrix_read_matrix_market (lunou, mat, gr, symm, sign)
    implicit none
    integer(ip)     , intent(in)            :: lunou
    type(fem_matrix), intent(inout)         :: mat
    type(fem_graph) , target, intent(inout) :: gr 
    integer(ip)     , intent(in), optional  :: symm
    integer(ip)     , intent(in), optional  :: sign

    integer(ip) :: i, j, iw1(2),iw2(2)
    integer(ip) :: nv, nw, nz

    integer(ip), allocatable :: ija_work(:,:), ija_index(:)
    real(rp)   , allocatable :: a_work(:)

    mat%storage =  scal    ! Scalar storage layout (blk: block; scal: scalar)
    mat%type    =  csr_mat ! fem_matrix type (csr_mat, csc_mat, css_mat)
    mat%nd1     =  1       ! Number of degrees of freedom, ndof1
    mat%nd2     =  1       ! Number of degrees of freedom, ndof2

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
          read(lunou,'(i10, i10, e32.25)',end=10, err=10) ija_work(1,i), ija_work(2,i), a_work(i)
          ija_index(i)=i
       end do
    else
       gr%type = csr_symm
       j=1
       do i=1,nz
          read(lunou,'(i10, i10, e32.25)', end=10,err=10) ija_work(1,j),ija_work(2,j),a_work(j)
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
    call memalloc (1,1,nz,mat%a,__FILE__,__LINE__)
    do i=1,nz
       mat%a(1,1,i)=a_work(ija_index(i))
    end do
    call memfree ( ija_index,__FILE__,__LINE__)
    call memfree (    a_work,__FILE__,__LINE__)

    mat%gr => gr

    return

10  write (0,*) 'Error reading matrix eof or err'
    stop

  end subroutine fem_matrix_read_matrix_market

  !=============================================================================
  subroutine fem_matrix_free_one_shot (f_matrix)
    implicit none
    type(fem_matrix), intent(inout) :: f_matrix

    call fem_matrix_free_progressively ( f_matrix, free_only_values )
    call fem_matrix_free_progressively ( f_matrix, free_only_struct )
    call fem_matrix_free_progressively ( f_matrix, free_clean )

  end subroutine fem_matrix_free_one_shot

    !=============================================================================
  subroutine fem_matrix_free_progressively (f_matrix, mode)
    implicit none
    type(fem_matrix), intent(inout) :: f_matrix
    integer(ip)     , intent(in)    :: mode 

    if ( mode == free_clean ) then
       f_matrix%nd1     = 0          ! Number of degrees of freedom
       f_matrix%nd2     = 0          ! Number of degrees of freedom
       f_matrix%storage = undef_sto

       ! AFM: bcs and msh are not required anymore
       ! in the new version of the codes
       ! nullify ( f_matrix%bcs )
       ! nullify ( f_matrix%msh )
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
       if(f_matrix%type==css_mat) then
          call memfree( f_matrix%d,__FILE__,__LINE__)
          call memfree( f_matrix%l,__FILE__,__LINE__)
          call memfree( f_matrix%u,__FILE__,__LINE__)
       else if(f_matrix%type==csr_mat) then
          call memfree( f_matrix%a,__FILE__,__LINE__)
       else if(f_matrix%type==csc_mat) then
          call memfree( f_matrix%a,__FILE__,__LINE__)
       end if
    end if

  end subroutine fem_matrix_free_progressively

  !=============================================================================
  subroutine fem_matrix_assembly_standard(nn,ln,ea,mat)
    implicit none
    integer(ip) , intent(in)    :: nn
    integer(ip) , intent(in)    :: ln(:)
    type(elmat) , intent(in)    :: ea
    type(fem_matrix), intent(inout) :: mat
    integer(ip) :: nd1,nd2,nl,ne

    nl=size(ln)
    assert(nl==nn)

    nd1=size(ea%data,dim=1)
    nd2=size(ea%data,dim=2)

    ! write (*,*) nd1, mat%nd1

    assert(nd1==mat%nd1)
    assert(nd2==mat%nd2)

    ne=size(ea%data,dim=3)
    assert(ne==nn)

    ne=size(ea%data,dim=4)
    assert(ne==nn)

    if (mat%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(mat%type==css_mat) then
          call ass_css_mat(mat%symm,nn,ln,ea%data,mat%nd1,mat%nd2,mat%gr%nv, &
               &           mat%gr%ia,mat%gr%is,mat%gr%ja,mat%d,mat%l,mat%u)
       else if(mat%type==csr_mat) then
          call ass_csr_mat(mat%symm,nn,ln,ea%data,mat%nd1,mat%nd2,mat%gr%nv, &
               &           mat%gr%ia,mat%gr%ja,mat%a)
       else if(mat%type==csc_mat) then
          call ass_csc_mat(mat%symm,nn,ln,ea%data,mat%nd1,mat%nd2,mat%gr%nv, &
               &           mat%gr%ia,mat%gr%ja,mat%a)
       end if
    else if (mat%storage == scal) then ! gr must describe scalar sparsity pattern of mat
       if( mat%type==css_mat ) then
          call ass_css_mat_scal(mat%symm,nn,ln,ea%data,mat%nd1,mat%nd2,mat%gr%nv, &
               &                mat%gr%ia,mat%gr%is,mat%gr%ja,mat%d,mat%l,mat%u)
       else if(mat%type==csr_mat) then
          call ass_csr_mat_scal(mat%symm,nn,ln,ea%data,mat%nd1,mat%nd2,mat%gr%nv, &
               &                mat%gr%ia,mat%gr%ja,mat%a)
       else if(mat%type==csc_mat) then
          call ass_csc_mat_scal(mat%symm,nn,ln,ea%data,mat%nd1,mat%nd2,mat%gr%nv, &
               &                mat%gr%ia,mat%gr%ja,mat%a)
       end if
    end if

  end subroutine fem_matrix_assembly_standard

  !=============================================================================
  subroutine fem_matrix_assembly_w_dof_coupl(nnode,nvars,co_dof_coupl,lst_coarse_dofs,lst_vars,elm,mat)
    implicit none
    integer(ip), intent(in)         :: nnode 
    integer(ip), intent(in)         :: nvars
    integer(ip), intent(in)         :: co_dof_coupl(nvars, nvars)
    integer(ip), intent(in)         :: lst_coarse_dofs(nnode)
    integer(ip), intent(in)         :: lst_vars(nnode) 
    type(elmat), intent(in)         :: elm 
    type(fem_matrix), intent(inout) :: mat 
    ! Locals
 
    if (mat%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(mat%type==css_mat) then

       else if(mat%type==csr_mat) then

       else if(mat%type==csc_mat) then

       end if
    else if (mat%storage == scal) then ! gr must describe scalar sparsity pattern of mat
       if( mat%type==css_mat ) then

       else if(mat%type==csr_mat) then
          call ass_csr_mat_scal_w_dof_handler(mat%symm, nnode, lst_coarse_dofs, lst_vars, elm%data, mat%gr%nv, &
               &                mat%gr%ia, mat%gr%ja, mat%a, nvars, co_dof_coupl)
          !call ass_csr_mat_scal_w_dof_handler(mat%symm, nn*dofh%nvarsxprob(el%prob), l2g1, ea%data, mat%gr%nv, &
          !     &                mat%gr%ia, mat%gr%ja, mat%a)
       else if(mat%type==csc_mat) then

       end if
    end if

  end subroutine fem_matrix_assembly_w_dof_coupl 

  !=============================================================================
  subroutine fem_matrix_assembly_w_dof_handler(nn,dofh,el,l2g1,l2g2,l2g3,l2g4,ea,mat)
    implicit none
    type(dof_handler), intent(in)    :: dofh
    type(fem_element), intent(in)    :: el
    type(elmat) ,      intent(in)    :: ea
    integer(ip) ,      intent(in)    :: nn(:)
    type(fem_matrix),  intent(inout) :: mat
    integer(ip),       intent(inout) :: l2g1(:),l2g2(:),l2g3(:),l2g4(:)
    ! Locals
    integer(ip) :: nd1,nd2,inode,i,j

    ! Asserts
    nd1=size(ea%data,dim=1)
    nd2=size(ea%data,dim=2)

    assert(nd1==mat%nd1)
    assert(nd2==mat%nd2)
    assert(nd1 == 1)
    assert(nd2 == 1)
    assert(size(nn,1)==el%nint) 

    l2g1 = 0; l2g2 = 0
    do j = dofh%i_varsxprob(1)%a(el%prob),dofh%i_varsxprob(1)%a(el%prob+1)-1
       do inode = 1,nn(el%iv(dofh%j_varsxprob(1)%a(j)))
          i = el%p_nod(inode) + j - dofh%i_varsxprob(1)%a(el%prob)
          l2g1(i) = el%elem2dof(inode,dofh%j_varsxprob(1)%a(j))  ! idof
          if(l2g1(i)/=0) l2g2(i) = dofh%j_varsxprob(1)%a(j)
       end do
    end do

    if (mat%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(mat%type==css_mat) then

       else if(mat%type==csr_mat) then

       else if(mat%type==csc_mat) then

       end if
    else if (mat%storage == scal) then ! gr must describe scalar sparsity pattern of mat
       if( mat%type==css_mat ) then

       else if(mat%type==csr_mat) then
          call ass_csr_mat_scal_w_dof_handler(mat%symm, el%ldof, l2g1, l2g2, ea%data, mat%gr%nv, &
               &                mat%gr%ia, mat%gr%ja, mat%a, dofh%nvars, dofh%dof_coupl, dofh)
          !call ass_csr_mat_scal_w_dof_handler(mat%symm, nn*dofh%nvarsxprob(el%prob), l2g1, ea%data, mat%gr%nv, &
          !     &                mat%gr%ia, mat%gr%ja, mat%a)
       else if(mat%type==csc_mat) then

       end if
    end if

  end subroutine fem_matrix_assembly_w_dof_handler

  !=============================================================================
  subroutine fem_matrix_assembly_nosq_w_dof_handler(nn,dofh,el,l2g1,l2g2,l2g3,l2g4,bl,ea,mat)
    implicit none
    integer(ip) ,      intent(in)    :: nn(:)
    type(dof_handler), intent(in)    :: dofh
    type(fem_element), intent(in)    :: el
    type(fem_blocks),  intent(in)    :: bl
    type(elmat) ,      intent(in)    :: ea
    type(fem_matrix),  intent(inout) :: mat
    integer(ip),       intent(inout) :: l2g1(:),l2g2(:),l2g3(:),l2g4(:)
    ! Locals
    integer(ip) :: nd1,nd2,nl,ne,inode,i,j,idofs,jdofs,maxnode

    ! Asserts
    assert(bl%nb==2)
    nd1=size(ea%data,dim=1)
    nd2=size(ea%data,dim=2)
    assert(nd1==mat%nd1)
    assert(nd2==mat%nd2)
    assert(nd1 == 1)
    assert(nd2 == 1)
    assert(size(nn,1)==el%nint) 

    maxnode=size(el%jvars,1)

    ! ROW DOFS & VARS
    idofs = 0; l2g1 = 0; l2g3 = 0
    do j = bl%ib(1),bl%ib(2)-1
       do inode = 1,nn(el%iv(bl%jb(j)))
          idofs = idofs+1
          l2g1(idofs) = el%elem2dof(inode,bl%jb(j))
          if(l2g1(idofs)/=0) l2g3(idofs) = bl%jb(j)     ! ivar
       end do
    end do

    ! COLUMN DOFS & VARS
    jdofs = 0; l2g2 = 0; l2g4 = 0
    do j = bl%ib(2),bl%ib(3)-1
       do inode = 1,nn(el%iv(bl%jb(j)))
          jdofs = jdofs+1
          l2g2(jdofs) = el%elem2dof(inode,bl%jb(j))
          if(l2g2(jdofs)/=0) l2g4(jdofs) = bl%jb(j)     ! jvar
       end do
    end do

    if (mat%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(mat%type==css_mat) then

       else if(mat%type==csr_mat) then

       else if(mat%type==csc_mat) then

       end if
    else if (mat%storage == scal) then ! gr must describe scalar sparsity pattern of mat
       if( mat%type==css_mat ) then

       else if(mat%type==csr_mat) then
          call ass_csr_mat_scal_nosq_w_dof_handler(mat%symm,el%nint,nn,size(el%p_nod,1),idofs,jdofs, &
               &                                   el%ldof,bl%ib,bl%jb,size(el%iv,1),el%iv,el%p_nod, &
               &                                   l2g1,l2g2,l2g3,l2g4,ea%data,mat%gr%nv,mat%gr%ia,  &
               &                                   mat%gr%ja,mat%a,dofh%nvars,dofh%dof_coupl,maxnode,&
               &                                   el%jvars)
       else if(mat%type==csc_mat) then

       end if
    end if

  end subroutine fem_matrix_assembly_nosq_w_dof_handler

  !=============================================================================
  subroutine fem_face_matrix_assembly_w_dof_handler(dofh,elem1,elem2,face,mat,elm,elm_aux,llnod,    &
       &                                            l2g,l2g2,ndime)

    implicit none
    type(dof_handler)  , intent(in)    :: dofh         ! Dof handler 
    type(fem_element)  , intent(in)    :: elem1, elem2 ! Neighbouring elements
    type(fem_face)     , intent(in)    :: face         ! Face
    type(fem_matrix)   , intent(inout) :: mat          ! Linear system matrix
    type(elmat)        , intent(inout) :: elm          ! Elemental matrix
    type(elmat_pointer), intent(in)    :: elm_aux(2)   ! Elemental matrix
    integer(ip)        , intent(inout) :: llnod(:), l2g(:), l2g2(:)
    integer(ip)        , intent(in)    :: ndime        ! Dimension of the problem

    integer(ip)              :: nnode, nodel(2), nodfa(2), nvars(2), ntx_i(2), ielem
    integer(ip)              :: inode, ivars,jnode,jvars, i
    !------------------------------------------------------------------------------------------
    ! When integrating on the faces we consider the nodes of the pair of elements neighboring
    ! that face. The element matrix is then of size 2*nnode x 2*nnode (the nodes on the face
    ! are duplicated in the continuous case but this is not a problem since the contributions
    ! of the different elements are then added on the assembling routine). This implies that
    ! some pairs of nodes in the elemental matrix may not be connected on the graph (for these
    ! elements the contribution of the elemental matrix is 0) thus the assembly can not be 
    ! performed directly but assembling element1 + face2 and element2+face1 separately.
    ! To do so we need to take into account that the terms of face1 + face2 are assembled
    ! twice so we must preprocess the matrix dividing those contributions by 2
    !------------------------------------------------------------------------------------------

    assert(elm%storage == scal)
    assert(elm_aux(1)%p%storage == elm%storage .and. elm_aux(1)%p%nd1 == elm%nd1 .and. elm_aux(1)%p%nd2 == elm%nd2)
    assert(elm_aux(2)%p%storage == elm%storage .and. elm_aux(2)%p%nd1 == elm%nd1 .and. elm_aux(2)%p%nd2 == elm%nd2)
    assert(elem1%prob==elem2%prob)
    assert(elem1%nint==1)

    !Set parameters:
    nodel(1) = elem1%f_inf(1)%p%nnode
    nodel(2) = elem2%f_inf(1)%p%nnode
    ntx_i(1) = elem1%f_inf(1)%p%ntxob_i(elem1%f_inf(1)%p%nobje_dim(ndime)+face%pos_elem(1)-1)
    ntx_i(2) = elem2%f_inf(1)%p%ntxob_i(elem2%f_inf(1)%p%nobje_dim(ndime)+face%pos_elem(2)-1)
    nodfa(1) = elem1%f_inf(1)%p%ntxob_i(elem1%f_inf(1)%p%nobje_dim(ndime)+face%pos_elem(1)) - ntx_i(1)
    nodfa(2) = elem2%f_inf(1)%p%ntxob_i(elem2%f_inf(1)%p%nobje_dim(ndime)+face%pos_elem(2)) - ntx_i(2)     
    nvars(1) = dofh%nvarsxprob(elem1%prob)
    nvars(2) = dofh%nvarsxprob(elem2%prob)


    !Treatment of the matrix so the terms of the face are not recomputed
    nnode = sum(nodfa)
    llnod(1:nodfa(1)) = elem1%f_inf(1)%p%ntxob_j(ntx_i(1):ntx_i(1)+nodfa(1)-1)
    llnod(nodfa(1)+1:nnode) = nodel(1) + elem2%f_inf(1)%p%ntxob_j(ntx_i(2):ntx_i(2)+nodfa(2)-1)
    do inode=1,nnode
       do ivars=1,elm%nd1
          do jnode=1,nnode
             do jvars=1,elm%nd2
                elm%data(1,1,(llnod(inode)-1)*elm%nd1+ivars,(llnod(jnode)-1)*elm%nd2+jvars) = &
                     & 0.5_rp * elm%data(1,1,(llnod(inode)-1)*elm%nd1+ivars,(llnod(jnode)-1)*elm%nd2+jvars)
             end do
          end do
       end do
    end do
    llnod = 0

    !Assembly matrix in two steps
    !FIRST STEP
    ielem = 1
    nnode = nodel(ielem) + nodfa(3-ielem)
    !Take the elm_aux submatrix to assemble
    llnod(1:nodel(ielem))       = (/(inode, inode=1,nodel(ielem))/) 
    llnod(nodel(ielem)+1:nnode) = nodel(ielem)+elem2%f_inf(1)%p%ntxob_j(ntx_i(2):ntx_i(2)+nodfa(2)-1)
    do inode=1,nnode
       do ivars=1,elm%nd1
          do jnode=1,nnode
             do jvars=1,elm%nd2
                elm_aux(ielem)%p%data(1,1,(inode-1)*elm_aux(ielem)%p%nd1+ivars, &
                     & (jnode-1)*elm_aux(ielem)%p%nd2+jvars) = &
                     &  elm%data(1,1,(llnod(inode)-1)*elm%nd1+ivars,(llnod(jnode)-1)*elm%nd2+jvars)
             end do
          end do
       end do
    end do

    !Local to global llnodl
    i = 0; l2g = 0; l2g2 = 0
    do inode = 1,nodel(ielem)
       do ivars = 1,nvars(ielem)
          i = i+1
          l2g(i) = elem1%elem2dof(inode,ivars)                                  ! idof
          if(l2g(i)/=0) l2g2(i) = ivars                                         ! ivar
       end do
    end do
    do inode = 1,nodfa(3-ielem)
       do ivars = 1,nvars(3-ielem)
          i = i+1
          l2g(i) = elem2%elem2dof(elem2%f_inf(1)%p%ntxob_j(ntx_i(2)+inode-1),ivars) !idof
          if(l2g(i)/=0) l2g2(i) = ivars                                         ! ivar
       end do
    end do

    nnode = nodel(ielem)*nvars(ielem) + nodfa(3-ielem)*nvars(3-ielem)
    if (mat%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(mat%type==css_mat) then

       else if(mat%type==csr_mat) then

       else if(mat%type==csc_mat) then

       end if
    else if (mat%storage == scal) then ! gr must describe scalar sparsity pattern of mat
       if( mat%type==css_mat ) then

       else if(mat%type==csr_mat) then
          call  ass_csr_mat_scal_w_dof_handler(mat%symm, nnode,l2g(1:nnode),l2g2(1:nnode),elm_aux(ielem)%p%data,&
               &                mat%gr%nv,mat%gr%ia, mat%gr%ja, mat%a, dofh%nvars, dofh%dof_coupl,dofh)
       else if(mat%type==csc_mat) then

       end if
    end if

    !SECOND STEP
    llnod = 0; l2g = 0; l2g2 = 0
    ielem = 2
    nnode = nodel(ielem) + nodfa(3-ielem)
    !Take the elm_aux submatrix to assemble
    llnod(1:nodel(ielem))       = nodel(3-ielem) + (/(inode, inode=1,nodel(ielem))/) 
    llnod(nodel(ielem)+1:nnode) = elem1%f_inf(1)%p%ntxob_j(ntx_i(1):ntx_i(1)+nodfa(1)-1)
    do inode=1,nnode
       do ivars=1,elm%nd1
          do jnode=1,nnode
             do jvars=1,elm%nd2
                elm_aux(ielem)%p%data(1,1,(inode-1)*elm_aux(ielem)%p%nd1+ivars, &
                     &                (jnode-1)*elm_aux(ielem)%p%nd2+jvars)   = &
                     &  elm%data(1,1,(llnod(inode)-1)*elm%nd1+ivars,(llnod(jnode)-1)*elm%nd2+jvars)
             end do
          end do
       end do
    end do

    !Local to global llnodl
    i = 0
    do inode = 1,nodel(ielem)
       do ivars = 1,nvars(ielem)
          i = i+1
          l2g(i) = elem2%elem2dof(inode,ivars)
          if(l2g(i)/=0) l2g2(i) = ivars                                         ! ivar
       end do
    end do
    do inode = 1,nodfa(3-ielem)
       do ivars = 1,nvars(3-ielem)
          i = i+1
          l2g(i) = elem1%elem2dof(elem1%f_inf(1)%p%ntxob_j(ntx_i(1)+inode-1),ivars)
          if(l2g(i)/=0) l2g2(i) = ivars                                         ! ivar
       end do
    end do
    nnode = nodel(ielem)*nvars(ielem) + nodfa(3-ielem)*nvars(3-ielem)
    if (mat%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(mat%type==css_mat) then

       else if(mat%type==csr_mat) then

       else if(mat%type==csc_mat) then

       end if
    else if (mat%storage == scal) then ! gr must describe scalar sparsity pattern of mat
       if( mat%type==css_mat ) then

       else if(mat%type==csr_mat) then
          call  ass_csr_mat_scal_w_dof_handler(mat%symm, nnode,l2g(1:nnode),l2g2(1:nnode),elm_aux(ielem)%p%data,&
               &                mat%gr%nv, mat%gr%ia, mat%gr%ja, mat%a, dofh%nvars, dofh%dof_coupl,dofh)
       else if(mat%type==csc_mat) then

       end if
    end if

  end subroutine fem_face_matrix_assembly_w_dof_handler

  !============================================================================

   subroutine add_contribution_to_global_matrix(n_glob_dofs, ia, ja, ig, jg, value_to_add, glob_matrix_a)

      implicit none
      integer(ip) , intent(in)      :: n_glob_dofs
      integer(ip) , intent(in)      :: ia(n_glob_dofs+1), ja(ia(n_glob_dofs+1)-1)
      integer(ip) , intent(in)      :: ig, jg       ! local and global indices
      real(rp), intent(in)          :: value_to_add
      real(rp)    , intent(inout)   :: glob_matrix_a(ia(n_glob_dofs+1)-1)

      integer(ip)                   :: k

      k = ia(ig)
      do while(ja(k) /= jg)
         assert(k < ia(ig + 1))
         k = k+1
      end do
      
      !write(*, "(A,I3,A,I3,A,F8.5,A,I3)") "adding to global matrix (", ig, ", ", jg, "), value ", value_to_add, ", k: ", k 
      glob_matrix_a(k) = glob_matrix_a(k) + value_to_add

   end subroutine add_contribution_to_global_matrix
  
  
   !============================================================================
   subroutine add_ced_contributions_to_global_matrix(n_glob_dofs, n_loc_dofs, block_i, block_j, ia, ja, il, jl, ig, jg, elem_loc_matrix_a, glob_matrix_a, dofhandler)

      implicit none
      integer(ip) , intent(in)      :: n_glob_dofs
      integer(ip) , intent(in)      :: n_loc_dofs
      integer(ip) , intent(in)      :: block_i, block_j
      integer(ip) , intent(in)      :: ia(n_glob_dofs+1), ja(ia(n_glob_dofs+1)-1)
      integer(ip) , intent(in)      :: il, jl, ig, jg       ! local and global indices
      real(rp)    , intent(in)      :: elem_loc_matrix_a(n_loc_dofs,n_loc_dofs)
      real(rp)    , intent(inout)   :: glob_matrix_a(ia(n_glob_dofs+1)-1)
      type(dof_handler), intent(in) :: dofhandler
      
      type(adapt_constraint), pointer :: constraint_i, constraint_j
      integer(ip)                     :: l_cing_i, l_cing_j, num_dofs_i, num_dofs_j, cing_glob_dof_i
      integer(ip)                     :: cing_glob_dof_j
      integer(ip)                     :: cing_node_i, cing_node_j, num_cing_nodes_i, num_cing_nodes_j
      real(rp)                        :: integ_value, coeff_i, coeff_j
      
      ! will it work otherwise? I do not know, has to be thought about
      assert(block_i == 1)
      assert(block_j == 1)
      
      if(ig < 0) then
         constraint_i => dofhandler%constraint_list%list(-ig)
         num_cing_nodes_i = constraint_i%num_cing_nodes
      else
         num_cing_nodes_i = 1      
      end if
      
      if(jg < 0) then
         assert(is_constrained(dofhandler%constraint_list, -jg, constraint_j))
         num_cing_nodes_j = constraint_j%num_cing_nodes
      else
         num_cing_nodes_j = 1      
      end if
      
      do l_cing_i = 1, num_cing_nodes_i       
         if(ig > 0) then
            cing_glob_dof_i = ig
            coeff_i = 1.0_rp
         else
            cing_node_i = constraint_i%cing_nodes(l_cing_i)
            coeff_i = constraint_i%cing_coefs(l_cing_i)
            num_dofs_i = dofhandler%ob2dof_p(block_i)%a(cing_node_i + 1) - dofhandler%ob2dof_p(block_i)%a(cing_node_i)
            assert(num_dofs_i <= 1)
            if(num_dofs_i == 0) then
               cycle
            end if
            ! todo: second column in ob2dof_l
            cing_glob_dof_i = dofhandler%ob2dof_l(block_i)%a(dofhandler%ob2dof_p(block_i)%a(cing_node_i), 1)
         end if
         
         do l_cing_j = 1, num_cing_nodes_j
            if(jg > 0) then
               cing_glob_dof_j = jg
               coeff_j = 1.0_rp
            else
               cing_node_j = constraint_j%cing_nodes(l_cing_j)
               coeff_j = constraint_j%cing_coefs(l_cing_j)
               num_dofs_j = dofhandler%ob2dof_p(block_j)%a(cing_node_j + 1) - dofhandler%ob2dof_p(block_j)%a(cing_node_j)
               assert(num_dofs_j <= 1)
               if(num_dofs_j == 0) then
                  cycle
               end if
               ! todo: second column in ob2dof_l
               cing_glob_dof_j = dofhandler%ob2dof_l(block_j)%a(dofhandler%ob2dof_p(block_j)%a(cing_node_j), 1)
            end if
            
            integ_value = elem_loc_matrix_a(il, jl) * coeff_i * coeff_j
            call add_contribution_to_global_matrix(n_glob_dofs, ia, ja, cing_glob_dof_i, cing_glob_dof_j, integ_value, glob_matrix_a)
         end do
      end do
      
   end subroutine add_ced_contributions_to_global_matrix   

  
  !============================================================================

  ! QUESTION (Pavel Kus) the following function is called (among others) from fem_matrix_assembly_w_dof_coupl, whrere 
  ! there is no dof_handler. Is it still used? If not, interface to the following may be simplified      
  subroutine ass_csr_mat_scal_w_dof_handler(symm, nl, l2g1, l2g2, ea, nv, ia, ja, la, nvar, dofc, dofhandler)

    implicit none
    integer(ip) , intent(in)      :: symm
    integer(ip) , intent(in)      :: nl, nv, nvar
    integer(ip) , intent(in)      :: l2g1(nl),l2g2(nl)
    integer(ip) , intent(in)      :: ia(nv+1), ja(ia(nv+1)-1)
    integer(ip) , intent(in)      :: dofc(nvar,nvar)
    real(rp)    , intent(in)      :: ea(nl,nl)
    real(rp)    , intent(inout)   :: la(ia(nv+1)-1)
    
    type(dof_handler), intent(in), optional   :: dofhandler

    ! local variables
    integer(ip)                   :: il, ig, jl, jg, k, iv, jv

    if ( symm == symm_false ) then
       ! CSR
       do il = 1,nl
          ig = l2g1(il)
          iv = l2g2(il)
          if (ig /= 0) then
             do jl = 1,nl
                jg = l2g1(jl)
                jv = l2g2(jl)
                if (jg /= 0) then
                   if(dofc(iv,jv)==1) then   
                      if((ig > 0) .and. (jg > 0)) then
                         !both are regular
                         call add_contribution_to_global_matrix(nv, ia, ja, ig, jg, ea(il, jl), la)
                      else
                         assert(present(dofhandler))
                         call add_ced_contributions_to_global_matrix(nv, nl, iv, jv, ia, ja, il, jl, ig, jg, ea, la, dofhandler)
                      end if
                   end if
                end if
             end do
          end if
       end do
       
    else
       ! CSR_SYMM
       do il = 1,nl
          ig = l2g1(il)
          iv = l2g2(il)
          if (ig /= 0) then
             do jl = 1,nl
                jg = l2g1(jl)
                jv = l2g2(jl)
                if (jg /= 0 .and. jg>=ig) then
                   if(dofc(iv,jv)==1) then
                      if((ig > 0) .and. (jg > 0)) then
                         !both are regular
                         call add_contribution_to_global_matrix(nv, ia, ja, ig, jg, ea(il, jl), la)
                      else
                         assert(present(dofhandler))
                         call add_ced_contributions_to_global_matrix(nv, nl, iv, jv, ia, ja, il, jl, ig, jg, ea, la, dofhandler)
                      end if
                   end if
                end if
             end do
          end if
       end do
       
    end if

    ! CSS
    ! idem for 
    ! jl=il,nl
    ! k = ia(min(ig,jg)) 
    ! do while(ja(k)/= max(ig,jg))

  end subroutine ass_csr_mat_scal_w_dof_handler

 !============================================================================
  subroutine ass_csr_mat_scal_nosq_w_dof_handler(symm,nint,nn,nd,id,jd,ld,ib,jb,nva,iv,pn, &
       &                                       l2g1,l2g2,l2g3,l2g4,ea,nv,ia,ja,la,nvar,dofc,&
       &                                       mn,jbn)

    implicit none
    integer(ip) , intent(in)      :: symm
    integer(ip) , intent(in)      :: nint, nv, nd, nva, id, jd, ld, nvar, mn
    integer(ip) , intent(in)      :: nn(nint),pn(nd),iv(nva)
    integer(ip) , intent(in)      :: ib(3),jb(ib(3)-1)
    integer(ip) , intent(in)      :: l2g1(id), l2g2(jd), l2g3(id), l2g4(jd)
    integer(ip) , intent(in)      :: ia(nv+1), ja(ia(nv+1)-1)
    integer(ip) , intent(in)      :: dofc(nvar,nvar)
    integer(ip) , intent(in)      :: jbn(mn,nva)
    real(rp)    , intent(in)      :: ea(ld,ld)
    real(rp)    , intent(inout)   :: la(ia(nv+1)-1)
!!$    integer(ip) , intent(in)      :: nn, nvp, nv
!!$    integer(ip) , intent(in)      :: ib(3),jb(ib(3)-1)
!!$    integer(ip) , intent(in)      :: l2g1(nn*(ib(2)-ib(1))), l2g2(nn*(ib(3)-ib(2)))
!!$    integer(ip) , intent(in)      :: ia(nv+1), ja(ia(nv+1)-1)
!!$    real(rp)    , intent(in)      :: ea(nn*nvp,nn*nvp)
!!$    real(rp)    , intent(inout)   :: la(ia(nv+1)-1)

    ! local variables
    integer(ip)                   :: in, jn, ip, jp, il, ig, jl, jg, k, ie, je, iva, jva

    assert ( symm == symm_false )

    ! CSR
    il = 0
    do ip = ib(1),ib(2)-1
       do in = 1,nn(iv(jb(ip)))
          il = il + 1
          ig = l2g1(il)
          iva = l2g3(il)
          if (ig /= 0) then
             ie = pn(in)-1 + jbn(in,jb(ip))
             jl = 0
             do jp = ib(2),ib(3)-1
                do jn = 1,nn(iv(jb(jp)))
                   jl = jl + 1
                   jg = l2g2(jl)
                   jva = l2g4(jl)
                   if (jg /= 0) then
                      if(dofc(iva,jva)==1) then
                         je = pn(jn)-1 + jbn(jn,jb(jp))
                         do k = ia(ig),ia(ig+1)-1
                            if(ja(k) == jg) then
                               la(k) = la(k) + ea(ie,je)
                               exit
                            end if
                         end do
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do

    ! CSS
    ! idem for 
    ! jl=il,nl
    ! k = ia(min(ig,jg)) 
    ! do while(ja(k)/= max(ig,jg))

  end subroutine ass_csr_mat_scal_nosq_w_dof_handler

  !=============================================================================
  subroutine fem_matrix_info ( f_mat, me, np )
    implicit none

    ! Parameters 
    type(fem_matrix)     , intent(in)    :: f_mat
    integer              , intent(out)   :: me
    integer              , intent(out)   :: np

    me = 0
    np = 1 
  end subroutine fem_matrix_info

  !=============================================================================
  subroutine fem_matrix_zero(mat)
    implicit none
    type(fem_matrix), intent(inout)       :: mat

    if(mat%type==css_mat) then
       mat%d = 0.0_rp
       mat%l = 0.0_rp
       if(mat%symm==symm_false) then
          mat%u = 0.0_rp
       end if
    else if(mat%type==csr_mat.or.mat%type==csc_mat) then
       mat%a = 0.0_rp
    end if
  end subroutine fem_matrix_zero

  !=============================================================================
  ! Subroutine to compute C <- A + alpha*B
  subroutine fem_matrix_sum(A,B,C,alpha)
    implicit none
    type(fem_matrix)  , intent(in)     :: A,B
    real(rp), optional, intent(in)     :: alpha
    type(fem_matrix)  , intent(inout)  :: C

    assert ( A%storage == B%storage  .and. A%storage == C%storage )
    assert ( A%type == B%type .and. A%type == C%type )
    assert ( A%symm == B%symm .and. A%symm == C%symm )
    assert ( A%nd1 == B%nd1 .and. A%nd1 == C%nd1 ) 
    assert ( A%nd2 == B%nd2 .and. A%nd2 == C%nd2 ) 

    if(present(alpha)) then
       if(C%type==css_mat) then
          C%d = A%d + alpha*B%d
          C%l = A%l + alpha*B%l
          if(C%symm==symm_false) then
             C%u = A%u + alpha*B%u
          end if
       else if(C%type==csr_mat.or.C%type==csc_mat) then
          C%a = A%a + alpha*B%a
       end if
    else
       if(C%type==css_mat) then
          C%d = A%d + B%d
          C%l = A%l + B%l
          if(C%symm==symm_false) then
             C%u = A%u + B%u
          end if
       else if(C%type==csr_mat.or.C%type==csc_mat) then
          C%a = A%a + B%a
       end if
    end if

  end subroutine fem_matrix_sum

  subroutine fem_matrix_transpose(A, A_t)

    type(fem_matrix), intent(in)     :: A         ! Input matrix
    type(fem_matrix), intent(out)    :: A_t       ! Output matrix


    ! Locals 
    type(fem_graph) :: aux_graph
    integer :: k,i,j

    if ( A%gr%type == csr ) then
       call fem_matrix_alloc ( scal, csr_mat, symm_false,                  &
            1, 1,     &
            A%gr, A_t )
    else if (A%gr%type == csr_symm) then
       call fem_matrix_alloc ( scal, csr_mat, symm_true,                   &
            1, 1,     &
            A%gr, A_t )
    end if

    aux_graph = A%gr

    if (A%gr%type == csr_symm) then 
       A_t%a(1,1,:) = A%a(1,1,:)
    elseif (A%gr%type == csr ) then
       k = 0    
       do i = 1, A_t%gr%nv
          do j=1, (A_t%gr%ia(i+1) - A_t%gr%ia(i) )
             k = k+1
             A_t%a(1,1,k) = A%a(1,1,aux_graph%ia( A_t%gr%ja(k) ) )
             aux_graph%ia(A_t%gr%ja(k)) = aux_graph%ia(A_t%gr%ja(k)) + 1
          end do
       end do
    end if

     end subroutine fem_matrix_transpose

end module fem_matrix_class
