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
module csr_sparse_matrix_names

USE types_names
USE memor_names
USE vector_names
USE serial_scalar_array_names
USE sparse_matrix_utils_names
USE sparse_matrix_parameters_names
USE base_sparse_matrix_names, only : base_sparse_matrix_t, base_sparse_matrix_iterator_t, coo_sparse_matrix_t

implicit none

# include "debug.i90"

private

    type, extends(base_sparse_matrix_t) :: csr_sparse_matrix_t
    private
        character(len=3)                 :: format_name = csr_format    ! String format id
        integer(ip)                      :: nnz = 0                     ! Number of non zeros
        integer(ip), allocatable, public :: irp(:)                      ! Row pointers
        integer(ip), allocatable, public :: ja(:)                       ! Column indices        
        real(rp),    allocatable, public :: val(:)                      ! Values
    contains
    private
        procedure, public :: is_by_rows                              => csr_sparse_matrix_is_by_rows
        procedure, public :: is_by_cols                              => csr_sparse_matrix_is_by_cols
        procedure, public :: get_format_name                         => csr_sparse_matrix_get_format_name
        procedure, public :: set_nnz                                 => csr_sparse_matrix_set_nnz
        procedure, public :: get_nnz                                 => csr_sparse_matrix_get_nnz
        procedure, public :: copy_to_coo_body                        => csr_sparse_matrix_copy_to_coo_body
        procedure, public :: copy_from_coo_body                      => csr_sparse_matrix_copy_from_coo_body
        procedure, public :: move_to_coo_body                        => csr_sparse_matrix_move_to_coo_body
        procedure, public :: move_from_coo_body                      => csr_sparse_matrix_move_from_coo_body
        procedure, public :: move_to_fmt_body                        => csr_sparse_matrix_move_to_fmt_body
        procedure, public :: move_from_fmt_body                      => csr_sparse_matrix_move_from_fmt_body
        procedure         :: allocate_numeric                        => csr_sparse_matrix_allocate_numeric
        procedure         :: allocate_symbolic                       => csr_sparse_matrix_allocate_symbolic
        procedure, public :: allocate_values_body                    => csr_sparse_matrix_allocate_values_body
        procedure, public :: initialize_values                       => csr_sparse_matrix_initialize_values
        procedure, public :: update_bounded_values_body              => csr_sparse_matrix_update_bounded_values_body
        procedure, public :: update_bounded_value_body               => csr_sparse_matrix_update_bounded_value_body
        procedure, public :: update_bounded_values_by_row_body       => csr_sparse_matrix_update_bounded_values_by_row_body
        procedure, public :: update_bounded_values_by_col_body       => csr_sparse_matrix_update_bounded_values_by_col_body
        procedure, public :: update_bounded_dense_values_body        => csr_sparse_matrix_update_bounded_dense_values_body
        procedure, public :: update_bounded_square_dense_values_body => csr_sparse_matrix_update_bounded_square_dense_values_body
        procedure, public :: update_values_body                      => csr_sparse_matrix_update_values_body
        procedure, public :: update_dense_values_body                => csr_sparse_matrix_update_dense_values_body
        procedure, public :: update_square_dense_values_body         => csr_sparse_matrix_update_square_dense_values_body
        procedure, public :: update_value_body                       => csr_sparse_matrix_update_value_body
        procedure, public :: update_values_by_row_body               => csr_sparse_matrix_update_values_by_row_body
        procedure, public :: update_values_by_col_body               => csr_sparse_matrix_update_values_by_col_body
        procedure, public :: split_2x2_symbolic                      => csr_sparse_matrix_split_2x2_symbolic
        procedure, public :: split_2x2_numeric                       => csr_sparse_matrix_split_2x2_numeric
        procedure         :: split_2x2_symbolic_body                 => csr_sparse_matrix_split_2x2_symbolic_body
        procedure         :: split_2x2_numeric_body                  => csr_sparse_matrix_split_2x2_numeric_body
        procedure, public :: permute_and_split_2x2_numeric           => csr_sparse_matrix_permute_and_split_2x2_numeric
        procedure, public :: permute_and_split_2x2_symbolic          => csr_sparse_matrix_permute_and_split_2x2_symbolic
        procedure         :: permute_and_split_2x2_symbolic_body     => csr_sparse_matrix_permute_and_split_2x2_symbolic_body
        procedure         :: permute_and_split_2x2_numeric_body      => csr_sparse_matrix_permute_and_split_2x2_numeric_body
        procedure         :: expand_matrix_numeric_array             => csr_sparse_matrix_expand_matrix_numeric_array
        procedure         :: expand_matrix_numeric_coo               => csr_sparse_matrix_expand_matrix_numeric_coo
        procedure         :: expand_matrix_symbolic_array            => csr_sparse_matrix_expand_matrix_symbolic_array
        procedure         :: expand_matrix_symbolic_coo              => csr_sparse_matrix_expand_matrix_symbolic_coo
        procedure         :: expand_matrix_numeric_array_body        => csr_sparse_matrix_expand_matrix_numeric_array_body
        procedure         :: expand_matrix_numeric_coo_body          => csr_sparse_matrix_expand_matrix_numeric_coo_body
        procedure         :: expand_matrix_symbolic_array_body       => csr_sparse_matrix_expand_matrix_symbolic_array_body
        procedure         :: expand_matrix_symbolic_coo_body         => csr_sparse_matrix_expand_matrix_symbolic_coo_body
        procedure, public :: extract_diagonal                        => csr_sparse_matrix_extract_diagonal
        procedure, public :: free_coords                             => csr_sparse_matrix_free_coords
        procedure, public :: free_val                                => csr_sparse_matrix_free_val
        procedure         :: apply_body                              => csr_sparse_matrix_apply_body
        procedure         :: apply_transpose_body                    => csr_sparse_matrix_apply_transpose_body
        procedure         :: apply_to_dense_matrix_body              => csr_sparse_matrix_apply_to_dense_matrix_body
        procedure         :: apply_transpose_to_dense_matrix_body    => csr_sparse_matrix_apply_transpose_to_dense_matrix_body
        procedure, public :: print_matrix_market_body                => csr_sparse_matrix_print_matrix_market_body
        procedure, public :: print                                   => csr_sparse_matrix_print
        procedure, public :: create_iterator                         => csr_sparse_matrix_create_iterator
        procedure, public :: get_entry                               => csr_sparse_matrix_get_entry
    end type csr_sparse_matrix_t

    !---------------------------------------------------------------------
    !< CSR MATRIX ITERATOR TYPE
    !---------------------------------------------------------------------
    type, extends(base_sparse_matrix_iterator_t) :: csr_sparse_matrix_iterator_t
       private
       integer(ip) :: row_index
       integer(ip) :: nnz_index
       type(csr_sparse_matrix_t), pointer :: matrix

     contains
       procedure, non_overridable :: create       => csr_sparse_matrix_iterator_create
       procedure                  :: init         => csr_sparse_matrix_iterator_init
       procedure                  :: free         => csr_sparse_matrix_iterator_free
       procedure                  :: next         => csr_sparse_matrix_iterator_next
       procedure                  :: has_finished => csr_sparse_matrix_iterator_has_finished
       procedure                  :: get_row      => csr_sparse_matrix_iterator_get_row
       procedure                  :: get_column   => csr_sparse_matrix_iterator_get_column
       procedure                  :: get_entry    => csr_sparse_matrix_iterator_get_entry
       procedure                  :: set_entry    => csr_sparse_matrix_iterator_set_entry
    end type csr_sparse_matrix_iterator_t

#ifdef ENABLE_MKL
    interface
        subroutine mkl_dcsrmm (transa, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc)
        import rp
        !-----------------------------------------------------------------
        ! http://software.intel.com/sites/products/documentation/hpc/
        ! compilerpro/en-us/cpp/win/mkl/refman/bla/functn_mkl_dcsrmm.html#functn_mkl_dcsrmm
        !-----------------------------------------------------------------
            character(len=1), intent(in)    :: transa
            integer,          intent(in)    :: m
            integer,          intent(in)    :: n
            integer,          intent(in)    :: k
            real(rp),         intent(in)    :: alpha
            character(len=*), intent(in)    :: matdescra
            real(rp),         intent(in)    :: val(*)
            integer,          intent(in)    :: indx(*)
            integer,          intent(in)    :: pntrb(m)
            integer,          intent(in)    :: pntre(m)
            real(rp),         intent(in)    :: b(ldb,*)
            integer,          intent(in)    :: ldb
            real(rp),         intent(in)    :: beta
            real(rp),         intent(inout) :: c(ldc,*)
            integer,          intent(in)    :: ldc
        end subroutine mkl_dcsrmm
    end interface
#endif

public :: csr_sparse_matrix_t
public :: csr_format

contains

    subroutine csr_sparse_matrix_set_nnz(this, nnz)
    !-----------------------------------------------------------------
    !< Get the number of non zeros of the matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nnz
    !-----------------------------------------------------------------
        this%nnz = nnz
    end subroutine csr_sparse_matrix_set_nnz


    function csr_sparse_matrix_get_nnz(this) result(nnz)
    !-----------------------------------------------------------------
    !< Get the number of non zeros of the matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in) :: this
        integer(ip)                            :: nnz
    !-----------------------------------------------------------------
        nnz = this%nnz
    end function csr_sparse_matrix_get_nnz


    function csr_sparse_matrix_is_by_rows(this) result(is_by_rows)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by rows
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in) :: this
        logical                                :: is_by_rows 
    !-----------------------------------------------------------------
        is_by_rows = .true.
    end function csr_sparse_matrix_is_by_rows


    function csr_sparse_matrix_is_by_cols(this) result(is_by_cols)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by columns
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in) :: this
        logical                                :: is_by_cols
    !-----------------------------------------------------------------
        is_by_cols = .false.
    end function csr_sparse_matrix_is_by_cols


    function csr_sparse_matrix_get_format_name(this) result(format_name)
    !-----------------------------------------------------------------
    !< Return a string with the sparse matrix format name
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in) :: this
        character(len=:), allocatable          :: format_name
    !-----------------------------------------------------------------
        format_name = this%format_name
    end function csr_sparse_matrix_get_format_name


    subroutine csr_sparse_matrix_copy_to_coo_body(this, to)
    !-----------------------------------------------------------------
    !< Copy this (CSR) -> to (COO)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in)    :: this
        type(coo_sparse_matrix_t),  intent(inout) :: to
        integer(ip)                               :: nnz
        integer(ip)                               :: i
        integer(ip)                               :: j
    !-----------------------------------------------------------------
        nnz = this%get_nnz()
        call to%free()
        if(this%get_num_rows() == this%get_num_cols()) then
            call to%create(num_rows_and_cols = this%get_num_rows(),          &
                           symmetric_storage = this%get_symmetric_storage(), &
                           is_symmetric      = this%is_symmetric(),          &
                           sign              = this%get_sign(),              &
                           nz                = nnz)
        else
            call to%create(num_rows = this%get_num_rows(), &
                           num_cols = this%get_num_cols() )
        endif
        call to%set_nnz(nnz)

        do i=1, this%get_num_rows()
            do j=this%irp(i),this%irp(i+1)-1
                to%ia(j)  = i
                to%ja(j)  = this%ja(j)
            end do
        end do
        if(.not. this%is_symbolic()) then
            call to%allocate_values_body(nnz)
            to%val(1:nnz) = this%val(1:nnz)
        endif
        call to%set_sort_status_by_rows()
        call to%set_state(this%get_state())
    end subroutine csr_sparse_matrix_copy_to_coo_body


    subroutine csr_sparse_matrix_copy_from_coo_body(this, from)
    !-----------------------------------------------------------------
    !< Copy from (COO) -> this (CSR)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        type(coo_sparse_matrix_t),  intent(in)    :: from
        type(coo_sparse_matrix_t)                 :: tmp
        integer(ip), allocatable                  :: itmp(:)
        integer(ip)                               :: nnz
        integer(ip)                               :: nr
        integer(ip)                               :: irw
        integer(ip)                               :: i
        integer(ip)                               :: j
    !-----------------------------------------------------------------
        call this%free()
        nr = from%get_num_rows()
        call this%set_num_rows(nr)
        call this%set_num_cols(from%get_num_cols())
        call this%set_symmetry(from%is_symmetric())
        call this%set_symmetric_storage(from%get_symmetric_storage())
        call this%set_sign(from%get_sign())
        if (.not. from%is_by_rows()) then 
            call tmp%copy_from_coo(from)
            call tmp%sort_and_compress()
            nnz = tmp%get_nnz()
            call this%set_nnz(nnz)
            call move_alloc(tmp%ia,itmp)
            call move_alloc(tmp%ja,this%ja)
            if(.not. tmp%is_symbolic()) call move_alloc(tmp%val,this%val)
            call tmp%free()
        else
            nnz = from%get_nnz()
            call this%set_nnz(nnz)
            call memalloc(nnz, itmp, __FILE__, __LINE__)
            call memalloc(nnz, this%ja, __FILE__, __LINE__)
            itmp            = from%ia(1:nnz)
            this%ja(1:nnz)  = from%ja(1:nnz)
            if(.not. from%is_symbolic()) then
                call memalloc(nnz, this%val, __FILE__, __LINE__)
                this%val(1:nnz) = from%val(1:nnz)
            endif
        endif
        call memalloc(this%get_num_cols()+1, this%irp, __FILE__, __LINE__)
        if(nnz <= 0) then
            this%irp(:) = 1
        else
            assert(nr>=itmp(nnz))
            this%irp(1) = 1

            j = 1 
            i = 1
            irw = itmp(j) ! sorted by rows

            outer: do 
                inner: do 
                    if (i >= irw) exit inner
                    assert(i<=nr) 
                    this%irp(i+1) = this%irp(i) 
                    i = i + 1
                end do inner

                j = j + 1
                if (j > nnz) exit
                if (itmp(j) /= irw) then 
                    this%irp(i+1) = j
                    irw = itmp(j) 
                    i = i + 1
                endif
                if (i>nr) exit
            enddo outer
            do 
                if (i>nr) exit
                this%irp(i+1) = j
                i = i + 1
            end do
        endif 
        call memfree(itmp, __FILE__, __LINE__)
        call this%set_state(from%get_state())
    end subroutine csr_sparse_matrix_copy_from_coo_body


    subroutine csr_sparse_matrix_move_to_coo_body(this, to)
    !-----------------------------------------------------------------
    !< Move this (CSR) -> to (COO)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        type(coo_sparse_matrix_t),  intent(inout) :: to
        integer(ip)                               :: nnz
        integer(ip)                               :: nr
        integer(ip)                               :: i
        integer(ip)                               :: j
    !-----------------------------------------------------------------
        nnz = this%get_nnz()
        nr = this%get_num_rows()
        call to%free()
        call to%set_num_rows(nr)
        call to%set_num_cols(this%get_num_cols())
        call to%set_symmetry(this%is_symmetric())
        call to%set_symmetric_storage(this%get_symmetric_storage())
        call to%set_sign(this%get_sign())
        call to%set_nnz(nnz)
        call memalloc(nnz, to%ia, __FILE__, __LINE__)
        call move_alloc(from=this%ja, to=to%ja)
        call move_alloc(from=this%val, to=to%val)
        do i=1, nr
            do j=this%irp(i),this%irp(i+1)-1
                to%ia(j)  = i
            end do
        end do
        call memfree(this%irp, __FILE__, __LINE__)
        call to%set_sort_status_by_rows()
        call to%set_state(this%get_state())
        call this%free()
    end subroutine csr_sparse_matrix_move_to_coo_body


    subroutine csr_sparse_matrix_move_from_coo_body(this, from)
    !-----------------------------------------------------------------
    !< Move from (COO) -> this (CSR)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        type(coo_sparse_matrix_t),  intent(inout) :: from
        integer(ip), allocatable                  :: itmp(:)
        integer(ip)                               :: nnz
        integer(ip)                               :: nr
        integer(ip)                               :: irw
        integer(ip)                               :: i
        integer(ip)                               :: j
    !-----------------------------------------------------------------
        call this%free()
        nr = from%get_num_rows()
        call this%set_num_rows(nr)
        call this%set_num_cols(from%get_num_cols())
        call this%set_symmetry(from%is_symmetric())
        call this%set_symmetric_storage(from%get_symmetric_storage())
        call this%set_sign(from%get_sign())
        if (.not. from%is_by_rows()) call from%sort_and_compress()
        nnz = from%get_nnz()
        call this%set_nnz(nnz)
        call move_alloc(from%ja,this%ja)
        call move_alloc(from%val,this%val)
        call memalloc(this%get_num_rows()+1, this%irp, __FILE__, __LINE__)
        if(nnz <= 0) then
            this%irp(:) = 1     
        else
            call move_alloc(from%ia,itmp)
            assert(nr>=itmp(nnz))
            this%irp(1) = 1

            j = 1 
            i = 1
            irw = itmp(j) ! sorted by rows

            outer: do 
                inner: do 
                    if (i >= irw) exit inner
                    assert(i<=nr) 
                    this%irp(i+1) = this%irp(i) 
                    i = i + 1
                end do inner

                j = j + 1
                if (j > nnz) exit
                if (itmp(j) /= irw) then 
                    this%irp(i+1) = j
                    irw = itmp(j) 
                    i = i + 1
                endif
                if (i>nr) exit
            enddo outer
            do 
                if (i>nr) exit
                this%irp(i+1) = j
                i = i + 1
            end do
            call memfree(itmp, __FILE__, __LINE__)
        endif 
        call this%set_state(from%get_state())
        call from%free()
    end subroutine csr_sparse_matrix_move_from_coo_body


    subroutine csr_sparse_matrix_move_to_fmt_body(this, to)
    !-----------------------------------------------------------------
    !< Move this (CRS) -> to (FMT)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),  intent(inout) :: this
        class(base_sparse_matrix_t), intent(inout) :: to
        type(coo_sparse_matrix_t)                  :: tmp
        integer(ip)                                :: nnz
        integer(ip)                                :: nr
    !-----------------------------------------------------------------
        select type (to)
            type is (coo_sparse_matrix_t) 
                call this%move_to_coo(to)
            type is (csr_sparse_matrix_t) 
                call to%free()
                call to%set_num_rows(this%get_num_cols())
                call to%set_num_cols(this%get_num_cols())
                call to%set_symmetry(this%is_symmetric())
                call to%set_symmetric_storage(this%get_symmetric_storage())
                call to%set_sign(this%get_sign())
                call to%set_nnz(this%get_nnz())
                call to%set_state(this%get_state())
                call move_alloc(this%irp, to%irp)
                call move_alloc(this%ja,  to%ja)
                call move_alloc(this%val, to%val)
                call this%free()
            class default
                call this%move_to_coo(tmp)
                call to%move_from_coo(tmp)
        end select
    end subroutine csr_sparse_matrix_move_to_fmt_body


    subroutine csr_sparse_matrix_move_from_fmt_body(this, from)
    !-----------------------------------------------------------------
    !< Move from (FMT) -> this (CSR)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),  intent(inout) :: this
        class(base_sparse_matrix_t), intent(inout) :: from
        type(coo_sparse_matrix_t)                  :: tmp
    !-----------------------------------------------------------------
        select type (from)
            type is (coo_sparse_matrix_t) 
                call this%move_from_coo(from)
            type is (csr_sparse_matrix_t)
                call this%set_num_rows(from%get_num_rows())
                call this%set_num_cols(from%get_num_cols())
                call this%set_symmetry(from%is_symmetric())
                call this%set_symmetric_storage(from%get_symmetric_storage())
                call this%set_sign(from%get_sign())
                call this%set_nnz(from%get_nnz())
                call this%set_state(from%get_state())
                call move_alloc(from%irp, this%irp)
                call move_alloc(from%ja,  this%ja)
                call move_alloc(from%val, this%val)
                call from%free()
            class default
                call from%move_to_coo(tmp)
                call this%move_from_coo(tmp)
        end select
    end subroutine csr_sparse_matrix_move_from_fmt_body


    subroutine matvec(num_rows, num_cols, irp, ja, val, alpha, x, beta, y)
    !-------------------------------------------------------------
    !< Sparse matrix vector product y = alpha*A*x + beta*y
    !-------------------------------------------------------------
        integer(ip), intent(in)    :: num_rows
        integer(ip), intent(in)    :: num_cols
        integer(ip), intent(in)    :: irp(num_rows+1)
        integer(ip), intent(in)    :: ja(irp(num_rows+1)-1)
        real(rp)   , intent(in)    :: val(irp(num_rows+1)-1)
        real(rp)   , intent(in)    :: alpha
        real(rp)   , intent(in)    :: x(num_cols)
        real(rp)   , intent(in)    :: beta
        real(rp)   , intent(inout) :: y(num_rows)
        integer(ip)                :: ir,ic, iz
    !-------------------------------------------------------------
        if(size(ja)==0) return
        do ir = 1, num_rows
           y(ir) = beta*y(ir)
           do iz = irp(ir), irp(ir+1)-1
              ic   = ja(iz)
              y(ir) = y(ir) + alpha*x(ic)*val(iz)
           end do ! iz
        end do ! ir
    end subroutine matvec


    subroutine transpose_matvec(num_rows, num_cols, irp, ja, val, alpha, x, beta, y)
    !-------------------------------------------------------------
    !< Transpose sparse matrix vector product y = alpha*A'*x + beta*y
    !-------------------------------------------------------------
        integer(ip), intent(in)    :: num_rows
        integer(ip), intent(in)    :: num_cols
        integer(ip), intent(in)    :: irp(num_rows+1)
        integer(ip), intent(in)    :: ja(irp(num_rows+1)-1)
        real(rp)   , intent(in)    :: val(irp(num_rows+1)-1)
        real(rp)   , intent(in)    :: alpha
        real(rp)   , intent(in)    :: x(num_rows)
        real(rp)   , intent(in)    :: beta
        real(rp)   , intent(inout) :: y(num_cols)
        integer(ip)                :: ir,ic, iz
    !-------------------------------------------------------------
        if(size(ja)==0) return
        y = beta*y
        do ir = 1, num_rows
           do iz = irp(ir), irp(ir+1)-1
              ic   = ja(iz)
              y(ic) = y(ic) + alpha*x(ir)*val(iz)
           end do ! iz
        end do ! ir
    end subroutine transpose_matvec


    subroutine matvec_symmetric_storage(num_rows, num_cols, irp, ja, val, alpha, x, beta, y)
    !-------------------------------------------------------------
    !< Symmetric stored sparse matrix vector product y = alpha*A*x + beta*y
    !-------------------------------------------------------------
        integer(ip), intent(in)    :: num_rows
        integer(ip), intent(in)    :: num_cols
        integer(ip), intent(in)    :: irp(num_rows+1)
        integer(ip), intent(in)    :: ja(irp(num_rows+1)-1)
        real(rp)   , intent(in)    :: val(irp(num_rows+1)-1)
        real(rp)   , intent(in)    :: alpha
        real(rp)   , intent(in)    :: x(num_cols)
        real(rp)   , intent(in)    :: beta
        real(rp)   , intent(inout) :: y(num_rows)
        integer(ip)                :: ir,ic, iz
    !-------------------------------------------------------------
        assert(num_rows==num_cols)
        if(size(ja)==0) return
        y = beta*y
        do ir = 1, num_rows
            y(ir) = y(ir) + x(ja(irp(ir)))*val(irp(ir))
            do iz = irp(ir)+1, irp(ir+1)-1
                ic = ja(iz)
                y(ir) = y(ir) + alpha*x(ic)*val(iz)
                y(ic) = y(ic) + alpha*x(ir)*val(iz)
            end do ! iz
        end do ! ir
    end subroutine matvec_symmetric_storage


    subroutine csr_sparse_matrix_apply_body(op,x,y) 
    !-----------------------------------------------------------------
    !< Apply matrix vector product y=op*x
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in)    :: op
        class(vector_t),            intent(in)    :: x
        class(vector_t) ,           intent(inout) :: y 
    !-----------------------------------------------------------------
        real(rp), pointer :: x_entries(:)
        real(rp), pointer :: y_entries(:)
        
        call x%GuardTemp()
        select type(x)
            class is (serial_scalar_array_t)
                select type(y)
                    class is(serial_scalar_array_t)
                        x_entries => x%get_entries()
                        y_entries => y%get_entries()
                        y_entries = 0.0_rp
                        if (op%get_symmetric_storage()) then
                            call matvec_symmetric_storage(            &
                                        num_rows = op%get_num_rows(), &
                                        num_cols = op%get_num_cols(), &
                                        irp      = op%irp,            &
                                        ja       = op%ja,             &
                                        val      = op%val,            &
                                        alpha    = 1.0_rp,            &
                                        x        = x_entries,         &
                                        beta     = 0.0_rp,            &
                                        y        = y_entries )
                        else
                            call matvec(num_rows = op%get_num_rows(), &
                                        num_cols = op%get_num_cols(), &
                                        irp      = op%irp,            &
                                        ja       = op%ja,             &
                                        val      = op%val,            &
                                        alpha    = 1.0_rp,            &
                                        x        = x_entries,         &
                                        beta     = 0.0_rp,            &
                                        y        = y_entries )
                    end if
                end select
            class DEFAULT
                check(.false.)
        end select
        call x%CleanTemp()
    end subroutine csr_sparse_matrix_apply_body


    subroutine csr_sparse_matrix_apply_transpose_body(op,x,y) 
    !-----------------------------------------------------------------
    !< Apply transpose matrix vector product y=op'*x
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in)    :: op
        class(vector_t),            intent(in)    :: x
        class(vector_t) ,           intent(inout) :: y 
    !-----------------------------------------------------------------
        real(rp), pointer :: x_entries(:)
        real(rp), pointer :: y_entries(:)
        
        call x%GuardTemp()
        select type(x)
            class is (serial_scalar_array_t)
                select type(y)
                    class is(serial_scalar_array_t)
                        x_entries => x%get_entries()
                        y_entries => y%get_entries()
                        y_entries = 0.0_rp
                        if (op%get_symmetric_storage()) then
                            call matvec_symmetric_storage(            &
                                        num_rows = op%get_num_rows(), &
                                        num_cols = op%get_num_cols(), &
                                        irp      = op%irp,            &
                                        ja       = op%ja,             &
                                        val      = op%val,            &
                                        alpha    = 1.0_rp,            &
                                        x        = x_entries,         &
                                        beta     = 0.0_rp,            &
                                        y        = y_entries )
                        else
                            call transpose_matvec(num_rows = op%get_num_rows(), &
                                        num_cols = op%get_num_cols(),           &
                                        irp      = op%irp,                      &
                                        ja       = op%ja,                       &
                                        val      = op%val,                      &
                                        alpha    = 1.0_rp,                      &
                                        x        = x_entries,                   &
                                        beta     = 0.0_rp,                      &
                                        y        = y_entries )
                    end if
                end select
            class DEFAULT
                check(.false.)
        end select
        call x%CleanTemp()
    end subroutine csr_sparse_matrix_apply_transpose_body


    subroutine csr_sparse_matrix_apply_to_dense_matrix_body(op, n, alpha, LDB, b, beta, LDC, c) 
    !-----------------------------------------------------------------
    !< Apply matrix matrix product y = alpha*op*b + beta*c
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),  intent(in)    :: op              ! Sparse matrix
        integer(ip),                 intent(in)    :: n               ! Number of columns of B and C dense arrays
        real(rp),                    intent(in)    :: alpha           ! Scalar alpha
        integer(ip),                 intent(in)    :: LDB             ! Leading dimensions of B matrix
        real(rp),                    intent(in)    :: b(LDB, n)       ! Matrix B
        real(rp),                    intent(in)    :: beta            ! Scalar beta
        integer(ip),                 intent(in)    :: LDC             ! Leading dimension of C matrix
        real(rp),                    intent(inout) :: c(LDC, n)       ! Matrix C
        integer(ip)                                :: i               ! Index to loop on B columns
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        if (op%get_symmetric_storage()) then
            call mkl_dcsrmm(transa    = 'N',               & ! Non transposed
                            m         = op%get_num_rows(), &
                            n         = n,                 &
                            k         = op%get_num_cols(), &
                            alpha     = alpha,             &
                            matdescra = 'SUNF',            & ! (Symmetric, Upper, Non-unit, Fortran)
                            val       = op%val,            &
                            indx      = op%ja,             &
                            pntrb     = op%irp(1),         &
                            pntre     = op%irp(2),         &
                            b         = b,                 &
                            ldb       = LDB,               &
                            beta      = beta,              &
                            c         = c,                 &
                            ldc       = LDC )
        else
            call mkl_dcsrmm(transa    = 'N',               & ! Non transposed
                            m         = op%get_num_rows(), &
                            n         = n,                 &
                            k         = op%get_num_cols(), &
                            alpha     = alpha,             &
                            matdescra = 'GXXF',            & ! General, X, X, Fortran)
                            val       = op%val,            &
                            indx      = op%ja,             &
                            pntrb     = op%irp(1),         &
                            pntre     = op%irp(2),         &
                            b         = b,                 &
                            ldb       = LDB,               &
                            beta      = beta,              &
                            c         = c,                 &
                            ldc       = LDC )
        endif
#else
        if (op%get_symmetric_storage()) then
            do i=1,n
                call matvec_symmetric_storage(                   &
                            num_rows = op%get_num_rows(),        &
                            num_cols = op%get_num_rows(),        &
                            irp      = op%irp,                   &
                            ja       = op%ja,                    &
                            val      = op%val,                   &
                            alpha    = alpha,                    &
                            x        = b(1:op%get_num_rows(),i), &
                            beta     = beta,                     &
                            y        = c(1:op%get_num_rows(),i) )
            end do
        else
            do i=1,n
                call matvec(num_rows = op%get_num_rows(),        &
                            num_cols = op%get_num_cols(),        &
                            irp      = op%irp,                   &
                            ja       = op%ja,                    &
                            val      = op%val,                   &
                            alpha    = alpha,                    &
                            x        = b(1:op%get_num_cols(),i), &
                            beta     = beta,                     &
                            y        = c(1:op%get_num_rows(),i) )
            enddo
        endif
#endif
    end subroutine csr_sparse_matrix_apply_to_dense_matrix_body


    subroutine csr_sparse_matrix_apply_transpose_to_dense_matrix_body(op, n, alpha, LDB, b, beta, LDC, c) 
    !-----------------------------------------------------------------
    !< Apply matrix matrix product y = alpha*op*b + beta*c
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),  intent(in)    :: op              ! Sparse matrix
        integer(ip),                 intent(in)    :: n               ! Number of columns of B and C dense arrays
        real(rp),                    intent(in)    :: alpha           ! Scalar alpha
        integer(ip),                 intent(in)    :: LDB             ! Leading dimensions of B matrix
        real(rp),                    intent(in)    :: b(LDB, n)       ! Matrix B
        real(rp),                    intent(in)    :: beta            ! Scalar beta
        integer(ip),                 intent(in)    :: LDC             ! Leading dimension of C matrix
        real(rp),                    intent(inout) :: c(LDC, n)       ! Matrix C
        integer(ip)                                :: i               ! Index to loop on B columns
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        if (op%get_symmetric_storage()) then
            call mkl_dcsrmm(transa    = 'T',               & ! Transposed
                            m         = op%get_num_rows(), &
                            n         = n,                 &
                            k         = op%get_num_cols(), &
                            alpha     = alpha,             &
                            matdescra = 'SUNF',            & ! (Symmetric, Upper, Non-unit, Fortran)
                            val       = op%val,            &
                            indx      = op%ja,             &
                            pntrb     = op%irp(1),         &
                            pntre     = op%irp(2),         &
                            b         = b,                 &
                            ldb       = LDB,               &
                            beta      = beta,              &
                            c         = c,                 &
                            ldc       = LDC )
        else
            call mkl_dcsrmm(transa    = 'T',               & ! Transposed
                            m         = op%get_num_rows(), &
                            n         = n,                 &
                            k         = op%get_num_cols(), &
                            alpha     = alpha,             &
                            matdescra = 'GXXF',            & ! General, X, X, Fortran)
                            val       = op%val,            &
                            indx      = op%ja,             &
                            pntrb     = op%irp(1),         &
                            pntre     = op%irp(2),         &
                            b         = b,                 &
                            ldb       = LDB,               &
                            beta      = beta,              &
                            c         = c,                 &
                            ldc       = LDC )
        endif
#else
        if (op%get_symmetric_storage()) then
            do i=1,n
                call matvec_symmetric_storage(                   &
                            num_rows = op%get_num_rows(),        &
                            num_cols = op%get_num_rows(),        &
                            irp      = op%irp,                   &
                            ja       = op%ja,                    &
                            val      = op%val,                   &
                            alpha    = alpha,                    &
                            x        = b(1:op%get_num_rows(),i), &
                            beta     = beta,                     &
                            y        = c(1:op%get_num_rows(),i) )
            end do
        else
            do i=1,n
                call transpose_matvec(                           &
                            num_rows = op%get_num_rows(),        &
                            num_cols = op%get_num_cols(),        &
                            irp      = op%irp,                   &
                            ja       = op%ja,                    &
                            val      = op%val,                   &
                            alpha    = alpha,                    &
                            x        = b(i,1:n), &
                            beta     = beta,                     &
                            y        = c(1:op%get_num_cols(),i) )
            enddo
        endif
#endif
    end subroutine csr_sparse_matrix_apply_transpose_to_dense_matrix_body


    subroutine csr_sparse_matrix_allocate_symbolic(this, nz)
    !-----------------------------------------------------------------
    !< Allocate coords
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip), optional,       intent(in)   :: nz
    !-----------------------------------------------------------------
        if(present(nz)) then
            call memalloc(nz, this%ja,  __FILE__, __LINE__)
        else
            call memalloc(max(7*this%get_num_rows(), 7*this%get_num_cols(), 1), this%ja,  __FILE__, __LINE__)
        endif
        call memalloc(this%get_num_rows()+1, this%irp,  __FILE__, __LINE__)
    end subroutine csr_sparse_matrix_allocate_symbolic


    subroutine csr_sparse_matrix_allocate_numeric(this, nz)
    !-----------------------------------------------------------------
    !< Allocate coords and values
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip), optional,       intent(in)   :: nz
    !-----------------------------------------------------------------
        call this%allocate_symbolic(nz)
        call this%allocate_values_body(nz)
    end subroutine csr_sparse_matrix_allocate_numeric


    subroutine csr_sparse_matrix_allocate_values_body(this, nz)
    !-----------------------------------------------------------------
    !< Allocate CSR values
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout)  :: this
        integer(ip), optional,      intent(in)     :: nz
    !-----------------------------------------------------------------
        assert(.not. allocated(this%val))
        if(present(nz)) then
            call memalloc(nz, this%val, __FILE__, __LINE__)
        else
            call memalloc(max(7*this%get_num_rows(), 7*this%get_num_cols(), 1), this%val,  __FILE__, __LINE__)
        endif
    end subroutine csr_sparse_matrix_allocate_values_body


    subroutine csr_sparse_matrix_initialize_values(this, val)
    !-----------------------------------------------------------------
    !< Initialize CSR values
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout)  :: this
        real(rp),                   intent(in)     :: val
    !-----------------------------------------------------------------
        if(allocated(this%val)) this%val(1:this%get_nnz()) = val
    end subroutine csr_sparse_matrix_initialize_values


    subroutine csr_sparse_matrix_update_bounded_values_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir,ic, ilr, ilc, ipaux,i1,i2,nr,nc
        logical                                   :: symmetric_storage
    !-----------------------------------------------------------------
        if(nz==0) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        ilr = -1 
        ilc = -1 
        symmetric_storage = this%get_symmetric_storage()
        do i=1, nz
            ir = ia(i)
            ic = ja(i) 
            ! Ignore out of bounds entries
            if (ir<imin .or. ir>imax .or. ic<jmin .or. ic>jmax .or. &
                ir<1 .or. ir>this%get_num_rows() .or. ic<1 .or. ic>this%get_num_cols() .or. &
                (symmetric_storage .and. ir>ic)) cycle
            i1 = this%irp(ir)
            i2 = this%irp(ir+1)
            nc = i2-i1
            ipaux = binary_search(ic,nc,this%ja(i1:i2-1))
            assert(ipaux>0) ! Entry not found
            if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
        end do
    end subroutine csr_sparse_matrix_update_bounded_values_body


    subroutine csr_sparse_matrix_update_bounded_value_body(this, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ipaux,i1,i2,nr,nc
    !-----------------------------------------------------------------
        ! Ignore out of bounds entries
        if (ia<imin .or. ia>imax .or. ja<jmin .or. ja>jmax .or. &
            ia<1 .or. ia>this%get_num_rows() .or. ja<1 .or. ja>this%get_num_cols() .or. &
            (this%get_symmetric_storage() .and. ia>ja)) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        i1 = this%irp(ia)
        i2 = this%irp(ia+1)
        nc = i2-i1
        ipaux = binary_search(ja,nc,this%ja(i1:i2-1))
        assert(ipaux>0) ! Entry not found
        if (ipaux>0) call apply_duplicates(input=val, output=this%val(i1+ipaux-1))
    end subroutine csr_sparse_matrix_update_bounded_value_body


    subroutine csr_sparse_matrix_update_bounded_values_by_row_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i, ic, ipaux,i1,i2,nc
        logical                                   :: symmetric_storage
    !-----------------------------------------------------------------
        if(nz==0 .or. ia<imin .or. ia<1 .or. ia>imax .or. ia>this%get_num_rows()) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        i1 = this%irp(ia)
        i2 = this%irp(ia+1)
        symmetric_storage = this%get_symmetric_storage()

        do i=1, nz
            ic = ja(i) 
            ! Ignore out of bounds entries
            if (ic<jmin .or. ic>jmax .or. ic<1 .or. ic>this%get_num_cols() .or. &
                (symmetric_storage .and. ia>ic)) cycle
            nc = i2-i1
            ipaux = binary_search(ic,nc,this%ja(i1:i2-1))
            assert(ipaux>0) ! Entry not found
            if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
        end do
    end subroutine csr_sparse_matrix_update_bounded_values_by_row_body



    subroutine csr_sparse_matrix_update_bounded_values_by_col_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir,ipaux,i1,i2,nc
        logical                                   :: symmetric_storage
    !-----------------------------------------------------------------
        if(nz==0 .or. ja<jmin .or. ja<1 .or. ja>jmax .or. ja>this%get_num_cols()) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        symmetric_storage = this%get_symmetric_storage()

        do i=1, nz
            ir = ia(i)
            ! Ignore out of bounds entries
            if (ir<imin .or. ir<1 .or. ir>imax .or. ir>this%get_num_rows() .or. &
                (symmetric_storage .and. ir>ja)) cycle
            i1 = this%irp(ir)
            i2 = this%irp(ir+1)
            nc = i2-i1
            ipaux = binary_search(ja,nc,this%ja(i1:i2-1))
            assert(ipaux>0) ! Entry not found
            if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
        end do
    end subroutine csr_sparse_matrix_update_bounded_values_by_col_body


    subroutine csr_sparse_matrix_update_values_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir,ic, ilr, ilc, ipaux,i1,i2,nr,nc
    !-----------------------------------------------------------------
        call this%update_body(nz, ia, ja , val, 1, this%get_num_rows(), 1, this%get_num_cols())
    end subroutine csr_sparse_matrix_update_values_body


    subroutine csr_sparse_matrix_update_bounded_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: num_cols
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_cols)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1 .or. num_cols<1) return
        do j=1, num_cols
            do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
            enddo
        enddo

    end subroutine csr_sparse_matrix_update_bounded_dense_values_body


    subroutine csr_sparse_matrix_update_bounded_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_rows)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1) return
        do j=1, num_rows
            do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
            enddo
        enddo
    end subroutine csr_sparse_matrix_update_bounded_square_dense_values_body


    subroutine csr_sparse_matrix_update_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: num_cols
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_cols)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1 .or. num_cols<1) return
        do j=1, num_cols  
           do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), 1, this%get_num_rows(), 1, this%get_num_cols())
            enddo
        enddo
    end subroutine csr_sparse_matrix_update_dense_values_body


    subroutine csr_sparse_matrix_update_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_rows)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1) return
        do j=1, num_rows
            do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), 1, this%get_num_rows(), 1, this%get_num_cols())
            enddo
        enddo
    end subroutine csr_sparse_matrix_update_square_dense_values_body


    subroutine csr_sparse_matrix_update_value_body(this, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ipaux,i1,i2,nr,nc
    !-----------------------------------------------------------------
        call this%update_body(ia, ja , val, 1, this%get_num_rows(), 1, this%get_num_cols())
    end subroutine csr_sparse_matrix_update_value_body


    subroutine csr_sparse_matrix_update_values_by_row_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ic, ipaux,i1,i2,nr,nc
    !-----------------------------------------------------------------
        call this%update_body(nz, ia, ja , val, 1, this%get_num_rows(), 1, this%get_num_cols())
    end subroutine csr_sparse_matrix_update_values_by_row_body


    subroutine csr_sparse_matrix_update_values_by_col_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val(nz)
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir,ic, ilr, ilc, ipaux,i1,i2,nr,nc
    !-----------------------------------------------------------------
        call this%update_body(nz, ia, ja , val, 1, this%get_num_rows(), 1, this%get_num_cols())
    end subroutine csr_sparse_matrix_update_values_by_col_body


    subroutine csr_sparse_matrix_split_2x2_numeric(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2
    !< A = [A_II A_IG]
    !<     [A_GI A_GG]
    !<
    !< this routine computes A_II, A_IG, A_GI and A_GG given the global 
    !< matrix A. Note that A_II, A_IG, A_GI and A_GG are all optional.
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        class(base_sparse_matrix_t),           intent(inout) :: A_II
        class(base_sparse_matrix_t),           intent(inout) :: A_IG
        class(base_sparse_matrix_t), optional, intent(inout) :: A_GI
        class(base_sparse_matrix_t),           intent(inout) :: A_GG
    !-----------------------------------------------------------------
        select type (A_II)
            type is (csr_sparse_matrix_t)
                select type(A_IG)
                    type is (csr_sparse_matrix_t)
                        select type(A_GG)
                            type is (csr_sparse_matrix_t)
                                if(present(A_GI)) then
                                    select type (A_GI)
                                        type is (csr_sparse_matrix_t)
                                            call this%split_2x2_numeric_body(num_row, num_col, A_II=A_II, A_IG=A_IG, A_GI=A_GI, A_GG=A_GG)
                                        class DEFAULT
                                            check(.false.)
                                    end select
                                else
                                    call this%split_2x2_numeric_body(num_row, num_col, A_II=A_II, A_IG=A_IG, A_GG=A_GG)
                                endif
                            class DEFAULT
                                check(.false.)
                        end select
                    class DEFAULT
                        check(.false.)
                end select
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_split_2x2_numeric


    subroutine csr_sparse_matrix_split_2x2_numeric_body(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 numeric implementation
    !< it must be called from split_2x2_symbolic where
    !< all preconditions are checked
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),          intent(in)    :: this
        integer(ip),                         intent(in)    :: num_row
        integer(ip),                         intent(in)    :: num_col
        type(csr_sparse_matrix_t),           intent(inout) :: A_II
        type(csr_sparse_matrix_t),           intent(inout) :: A_IG
        type(csr_sparse_matrix_t), optional, intent(inout) :: A_GI
        type(csr_sparse_matrix_t),           intent(inout) :: A_GG
        integer(ip)                                        :: i, j, k
        integer(ip)                                        :: nz
        integer(ip)                                        :: counter
        integer(ip)                                        :: row_index
        integer(ip)                                        :: nz_offset
        integer(ip)                                        :: nz_ignored
        integer(ip)                                        :: A_XX_lbound
        integer(ip)                                        :: A_XX_ubound
        integer(ip)                                        :: this_lbound
        integer(ip)                                        :: this_ubound
        integer(ip)                                        :: total_cols
        integer(ip)                                        :: total_rows
        integer(ip)                                        :: sign
        integer(ip)                                        :: state
        logical                                            :: is_properties_set_state
    !-----------------------------------------------------------------
        ! Check state
        assert(this%state_is_assembled()) 
        state = A_II%get_state() 
        is_properties_set_state = A_II%state_is_properties_set() 
        assert(state == A_IG%get_state() .and.  state == A_GG%get_state())
        if(present(A_GI)) then
            assert(state == A_II%get_state())
        endif
        assert(is_properties_set_state  .or. A_II%state_is_assembled_symbolic())

        total_rows = this%get_num_rows()
        total_cols = this%get_num_cols()

        if(is_properties_set_state) then
            ! Set properties to all submatrices
            call A_II%set_num_rows(num_row); call A_II%set_num_cols(num_col)
            call A_IG%set_num_rows(num_row); call A_IG%set_num_cols(total_cols-num_col)
            if(present(A_GI)) then
                assert(A_GI%get_state() == state)
                call A_GI%set_num_rows(total_rows-num_row); call A_GI%set_num_cols(num_col)
            endif
            call A_GG%set_num_rows(total_rows-num_row); call A_GG%set_num_cols(total_cols-num_col)
    
            ! Allocate irp, ja and val arrays of all submatrices
            nz = this%irp(num_row+1)-1  ! nnz after num_row
            if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
                call A_II%allocate_numeric(2*nz)
            else
                call A_II%allocate_numeric(nz)
            endif
            A_II%irp(1) = 1
            call A_IG%allocate_numeric(nz);     A_IG%irp(1) = 1
            nz = this%nnz - nz          ! nnz after num_row
            if(present(A_GI)) then
                call A_GI%allocate_numeric(nz); A_GI%irp(1) = 1
            endif
            if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
                call A_GG%allocate_numeric(2*nz)
            else
                call A_GG%allocate_numeric(nz)
            endif
            A_GG%irp(1) = 1
        else
            call A_II%allocate_values_body(A_II%nnz)
            call A_IG%allocate_values_body(A_IG%nnz)
            if(present(A_GI)) call A_GI%allocate_values_body(A_GI%nnz)
            call A_GG%allocate_values_body(A_GG%nnz)
        endif

        if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
            ! If input matrix has symmetric storage and one or more of the output diagonal matrix has not symmetric 
            A_II%irp = 0
            ! Count total number of
            do i=1, num_row
                do j=this%irp(i), this%irp(i+1)-1
                    k = this%ja(j)
                    if(k>num_col) exit
                    A_II%irp(i+1) = A_II%irp(i+1)+1
                    if(i == k) cycle ! diagonal elements
                    A_II%irp(k+1) = A_II%irp(k+1)+1
                enddo
            enddo
            ! Build CSR header from the counter
            A_II%irp(1) = 1
            do i=1, num_row
                A_II%irp(i+1) = A_II%irp(i+1)+A_II%irp(i)
            enddo
        endif

        ! Loop (1:num_row) to get A_II and A_IG 
        do i=1, num_row
            nz_offset = 0
            ! Count the number of columns less than or equal to num_row
            do j=this%irp(i), this%irp(i+1)-1
                if(this%ja(j)> num_col) exit
                nz_offset = nz_offset + 1
            enddo

            ! Number of nnz of A_II in row i
            nz = nz_offset 
            if(this%get_symmetric_storage() .eqv. A_II%get_symmetric_storage()) then
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_II%irp(i); A_XX_ubound = A_II%irp(i)+nz-1
                    this_lbound = this%irp(i); this_ubound = this%irp(i)+nz-1
                    ! Assign columns and values
                    if(is_properties_set_state) then
                        A_II%irp(i+1)                     = A_XX_ubound+1
                        A_II%ja (A_XX_lbound:A_XX_ubound) = this%ja(this_lbound:this_ubound)
                        A_II%nnz = A_II%nnz + nz
                    endif
                    A_II%val (A_XX_lbound:A_XX_ubound) = this%val(this_lbound:this_ubound)
                else
                    A_II%irp(i+1) = A_II%irp(i)
                endif

            else if(.not. this%get_symmetric_storage() .and. A_II%get_symmetric_storage()) then
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_II%irp(i); A_XX_ubound = A_II%irp(i)+nz-1
                    this_lbound = this%irp(i); this_ubound = this%irp(i)+nz-1
                    ! Ignore lower triangle entries
                    nz_ignored = 0
                    do j=this_lbound,this_ubound
                        if(this%ja(j)>=i) exit
                        nz_ignored = nz_ignored+1
                    enddo
                    ! Assign columns and values
                    if(is_properties_set_state) then
                        A_II%irp(i+1)                                = A_XX_ubound-nz_ignored+1
                        A_II%ja(A_XX_lbound:A_XX_ubound-nz_ignored) = this%ja(this_lbound+nz_ignored:this_ubound)
                        A_II%nnz = A_II%nnz - nz_ignored + nz
                    endif
                    A_II%val(A_XX_lbound:A_XX_ubound-nz_ignored) = this%val(this_lbound+nz_ignored:this_ubound)
                else
                    A_II%irp(i+1) = A_II%irp(i)
                endif

            else if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
                do j=this%irp(i), this%irp(i)+nz-1
                    ! Add upper triangle elements
                    k = this%ja(j)
                    A_II%ja(A_II%irp(i)) = k
                    A_II%val(A_II%irp(i)) = this%val(j)
                    A_II%irp(i) = A_II%irp(i)+1
                    if(i==k) cycle
                    ! Add lower triangle elements
                    A_II%ja(A_II%irp(k)) = i
                    A_II%val(A_II%irp(k)) = this%val(i)
                    A_II%irp(k) = A_II%irp(k)+1
                enddo
            endif

            ! Number of nnz of A_IG in row i
            nz = this%irp(i+1)-this%irp(i)-nz_offset 
            if(nz>0) then
                ! Calculate bounds
                A_XX_lbound = A_IG%irp(i);           A_XX_ubound = A_IG%irp(i)+nz-1
                this_lbound = this%irp(i)+nz_offset; this_ubound = this%irp(i+1)-1
                ! Assign columns and values
                if(is_properties_set_state) then
                    A_IG%irp(i+1)                     = A_XX_ubound+1
                    A_IG%ja (A_XX_lbound:A_XX_ubound) = this%ja (this_lbound:this_ubound)-num_col
                    A_IG%nnz = A_IG%nnz + nz
                endif
                A_IG%val(A_XX_lbound:A_XX_ubound) = this%val(this_lbound:this_ubound)
            else
                A_IG%irp(i+1) = A_IG%irp(i)
            endif
        enddo

        if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
            do i=num_row, 1, -1
                A_II%irp(i+1) = A_II%irp(i)
            enddo
            A_II%irp(1) = 1
            A_II%nnz = A_II%irp(num_row+1)-1
        endif

        if(is_properties_set_state) then
            call memrealloc(A_II%nnz, A_II%ja,   __FILE__, __LINE__)
            call memrealloc(A_IG%nnz, A_IG%ja,   __FILE__, __LINE__)
            call memrealloc(A_II%nnz, A_II%val,  __FILE__, __LINE__)
            call memrealloc(A_IG%nnz, A_IG%val,  __FILE__, __LINE__)
        endif

        if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
            ! If input matrix has symmetric storage and one or more of the output diagonal matrix has not symmetric 
            A_GG%irp = 0
            ! Count total number of
            do i=num_row+1, total_rows
                do j=this%irp(i), this%irp(i+1)-1
                    k = this%ja(j)-num_col
                    if(k<0) cycle
                    A_GG%irp(i-num_row+1) = A_GG%irp(i-num_row+1)+1
                    if(i-num_row == k) cycle ! diagonal elements
                    A_GG%irp(k+1) = A_GG%irp(k+1)+1
                enddo
            enddo
            ! Build CSR header from counter
            A_GG%irp(1) = 1
            do i=1, A_GG%get_num_rows()
                A_GG%irp(i+1) = A_GG%irp(i+1)+A_GG%irp(i)
            enddo
        endif

        ! Loop (num_row:this%num_rows) to get A_GI and A_GG
        do i=num_row+1, this%get_num_rows()
            nz_offset = 0
            ! Count the number of columns less than or equal to num_row
            do j=this%irp(i), this%irp(i+1)-1
                if(this%ja(j)> num_col) exit
                nz_offset = nz_offset + 1
            enddo

            ! Number of nnz of A_GI in row i-num_row
            nz = nz_offset 
            if(present(A_GI)) then
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_GI%irp(i-num_row); A_XX_ubound = A_GI%irp(i-num_row)+nz-1
                    this_lbound = this%irp(i);         this_ubound = this%irp(i)+nz-1
                    ! Assign columns and values
                    if(is_properties_set_state) then
                        A_GI%irp(i-num_row+1)             = A_XX_ubound+1
                        A_GI%ja (A_XX_lbound:A_XX_ubound) = this%ja (this_lbound:this_ubound)
                        A_GI%nnz = A_GI%nnz + nz
                    endif
                    A_GI%val(A_XX_lbound:A_XX_ubound) = this%val(this_lbound:this_ubound)
                else
                    A_GI%irp(i-num_row+1) = A_GI%irp(i-num_row)
                endif
            endif

            ! Number of nnz of A_GG in row i-num_row
            nz = this%irp(i+1)-this%irp(i)-nz_offset 
            if(this%get_symmetric_storage() .eqv. A_GG%get_symmetric_storage()) then
                ! Number of nnz of A_GG in row i-num_row
                nz = this%irp(i+1)-this%irp(i)-nz_offset 
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_GG%irp(i-num_row);   A_XX_ubound = A_GG%irp(i-num_row)+nz-1
                    this_lbound = this%irp(i)+nz_offset; this_ubound = this%irp(i+1)-1
                    ! Assign columns and values
                    if(is_properties_set_state) then
                        A_GG%irp(i-num_row+1)             = A_XX_ubound+1
                        A_GG%ja (A_XX_lbound:A_XX_ubound) = this%ja (this_lbound:this_ubound)-num_col
                        A_GG%nnz = A_GG%nnz + nz
                    endif
                    A_GG%val(A_XX_lbound:A_XX_ubound) = this%val(this_lbound:this_ubound)
                else
                    if(is_properties_set_state) A_GG%irp(i-num_row+1) = A_GG%irp(i-num_row)
                endif

            else if(.not. this%get_symmetric_storage() .and. A_GG%get_symmetric_storage()) then

                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_GG%irp(i-num_row);   A_XX_ubound = A_GG%irp(i-num_row)+nz-1
                    this_lbound = this%irp(i)+nz_offset; this_ubound = this%irp(i+1)-1
                    ! Ignore lower triangle entries
                    nz_ignored = 0
                    do j=this_lbound,this_ubound
                        if(this%ja(j)>=i) exit
                        nz_ignored = nz_ignored+1
                    enddo
                    ! Assign columns and values
                    if(is_properties_set_state) then
                        A_GG%irp(i-num_row+1)             = A_XX_ubound-nz_ignored+1
                        A_GG%ja (A_XX_lbound:A_XX_ubound-nz_ignored) = this%ja (this_lbound+nz_ignored:this_ubound)-num_col
                        A_GG%nnz = A_GG%nnz - nz_ignored + nz
                    endif
                    A_GG%val(A_XX_lbound:A_XX_ubound-nz_ignored) = this%val(this_lbound+nz_ignored:this_ubound)
                else
                    if(is_properties_set_state) A_GG%irp(i-num_row+1) = A_GG%irp(i-num_row)
                endif

            else if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
                do j=this%irp(i), this%irp(i)+nz-1
                    ! Add upper triangle elements
                    k = this%ja(j)-num_col
                    A_GG%ja(A_GG%irp(i-num_row)) = k
                    A_GG%val(A_GG%irp(i-num_row)) = this%val(j)
                    A_GG%irp(i-num_row) = A_GG%irp(i-num_row)+1
                    if(i-num_row==k) cycle
                    ! Add lower triangle elements
                    A_GG%ja(A_GG%irp(k)) = i-num_row
                    A_GG%val(A_GG%irp(k)) = this%val(i)
                    A_GG%irp(k) = A_GG%irp(k)+1
                enddo
            endif
        enddo

        if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
            do i=A_GG%get_num_rows(), 1, -1
                A_GG%irp(i+1) = A_GG%irp(i)
            enddo
            A_GG%irp(1) = 1
            A_GG%nnz = A_GG%irp(A_GG%get_num_rows()+1)-1
        endif

        if(present(A_GI) .and. is_properties_set_state) then
            call memrealloc(A_GI%nnz, A_GI%ja,   __FILE__, __LINE__)
            call memrealloc(A_GI%nnz, A_GI%val,  __FILE__, __LINE__)
        endif
        if(is_properties_set_state)  then
            call memrealloc(A_GG%nnz, A_GG%ja,   __FILE__, __LINE__)
            call memrealloc(A_GG%nnz, A_GG%val,  __FILE__, __LINE__)
        endif

        call A_II%set_state_assembled()
        call A_IG%set_state_assembled()
        if(present(A_GI)) call A_GI%set_state_assembled()
        call A_GG%set_state_assembled()
    end subroutine csr_sparse_matrix_split_2x2_numeric_body


    subroutine csr_sparse_matrix_split_2x2_symbolic(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2
    !< A = [A_II A_IG]
    !<     [A_GI A_GG]
    !<
    !< this routine computes A_II, A_IG, A_GI and A_GG given the global 
    !< matrix A. Note that A_II, A_IG, A_GI and A_GG are all optional.
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        class(base_sparse_matrix_t),           intent(inout) :: A_II
        class(base_sparse_matrix_t),           intent(inout) :: A_IG
        class(base_sparse_matrix_t), optional, intent(inout) :: A_GI
        class(base_sparse_matrix_t),           intent(inout) :: A_GG
    !-----------------------------------------------------------------
        select type (A_II)
            type is (csr_sparse_matrix_t)
                select type(A_IG)
                    type is (csr_sparse_matrix_t)
                        select type(A_GG)
                            type is (csr_sparse_matrix_t)
                                if(present(A_GI)) then
                                    select type (A_GI)
                                        type is (csr_sparse_matrix_t)
                                            call this%split_2x2_symbolic_body(num_row, num_col, A_II=A_II, A_IG=A_IG, A_GI=A_GI, A_GG=A_GG)
                                        class DEFAULT
                                            check(.false.)
                                    end select
                                else
                                    call this%split_2x2_symbolic_body(num_row, num_col, A_II=A_II, A_IG=A_IG, A_GG=A_GG)
                                endif
                            class DEFAULT
                                check(.false.)
                        end select
                    class DEFAULT
                        check(.false.)
                end select
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_split_2x2_symbolic


    subroutine csr_sparse_matrix_split_2x2_symbolic_body(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 symbolic implementation
    !< it must be called from split_2x2_symbolic where
    !< all preconditions are checked
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),          intent(in)    :: this
        integer(ip),                         intent(in)    :: num_row
        integer(ip),                         intent(in)    :: num_col
        type(csr_sparse_matrix_t),           intent(inout) :: A_II
        type(csr_sparse_matrix_t),           intent(inout) :: A_IG
        type(csr_sparse_matrix_t), optional, intent(inout) :: A_GI
        type(csr_sparse_matrix_t),           intent(inout) :: A_GG
        integer(ip)                                        :: i, j, k
        integer(ip)                                        :: nz
        integer(ip)                                        :: row_index
        integer(ip)                                        :: nz_offset
        integer(ip)                                        :: nz_ignored
        integer(ip)                                        :: A_XX_lbound
        integer(ip)                                        :: A_XX_ubound
        integer(ip)                                        :: this_lbound
        integer(ip)                                        :: this_ubound
        integer(ip)                                        :: total_cols
        integer(ip)                                        :: total_rows
        integer(ip)                                        :: sign
        integer(ip)                                        :: state
        logical                                            :: properties_are_setted
    !-----------------------------------------------------------------
        ! Check state
        assert(this%state_is_assembled() .or. this%state_is_assembled_symbolic() ) 
        state = A_II%get_state()
        properties_are_setted = A_II%state_is_properties_set()
        assert( state == A_IG%get_state() .and. state== A_GG%get_state())
        assert(properties_are_setted)

        total_rows = this%get_num_rows()
        total_cols = this%get_num_cols()
        assert(num_row>0 .and. num_row<total_rows .and. num_col>0 .and. num_col<total_cols)

        call A_II%set_num_rows(num_row); call A_II%set_num_cols(num_col)
        call A_IG%set_num_rows(num_row); call A_IG%set_num_cols(total_cols-num_col)
        if(present(A_GI)) then
            assert(A_GI%get_state() == state)
            call A_GI%set_num_rows(total_rows-num_row); call A_GI%set_num_cols(num_col)
        endif
        call A_GG%set_num_rows(total_rows-num_row); call A_GG%set_num_cols(total_cols-num_col)

        ! Allocate irp, ja and val arrays of all submatrices
        nz = this%irp(num_row+1)-1  ! nnz after num_row
        if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
            call A_II%allocate_symbolic(2*nz)
        else
            call A_II%allocate_symbolic(nz)
        endif
        A_II%irp(1) = 1
        call A_IG%allocate_symbolic(nz);     A_IG%irp(1) = 1
        nz = this%nnz - nz          ! nnz after num_row
        if(present(A_GI)) then
            call A_GI%allocate_symbolic(nz); A_GI%irp(1) = 1
        endif
        if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
            call A_GG%allocate_symbolic(2*nz)
        else
            call A_GG%allocate_symbolic(nz)
        endif
        A_GG%irp(1) = 1

        if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
            ! If input matrix has symmetric storage and one or more of the output diagonal matrix has not symmetric 
            A_II%irp = 0
            ! Count total number of
            do i=1, num_row
                do j=this%irp(i), this%irp(i+1)-1
                    k = this%ja(j)
                    if(k>num_col) exit
                    A_II%irp(i+1) = A_II%irp(i+1)+1
                    if(i == k) cycle ! diagonal elements
                    A_II%irp(k+1) = A_II%irp(k+1)+1
                enddo
            enddo
            A_II%irp(1) = 1
            do i=1, num_row
                A_II%irp(i+1) = A_II%irp(i+1)+A_II%irp(i)
            enddo
        endif

        ! Loop (1:num_row) to get A_II and A_IG 
        do i=1, num_row
            nz_offset = 0
            ! Count the number of columns less than or equal to num_row
            do j=this%irp(i), this%irp(i+1)-1
                if(this%ja(j)> num_col) exit
                nz_offset = nz_offset + 1
            enddo

            ! Number of nnz of A_II in row i
            nz = nz_offset 
            if(this%get_symmetric_storage() .eqv. A_II%get_symmetric_storage()) then
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_II%irp(i); A_XX_ubound = A_II%irp(i)+nz-1
                    this_lbound = this%irp(i); this_ubound = this%irp(i)+nz-1
                    ! Assign columns
                    A_II%irp(i+1)                     = A_XX_ubound+1
                    A_II%ja (A_XX_lbound:A_XX_ubound) = this%ja(this_lbound:this_ubound)
                    A_II%nnz = A_II%nnz + nz
                else
                    A_II%irp(i+1) = A_II%irp(i)
                endif

            else if(.not. this%get_symmetric_storage() .and. A_II%get_symmetric_storage()) then
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_II%irp(i); A_XX_ubound = A_II%irp(i)+nz-1
                    this_lbound = this%irp(i); this_ubound = this%irp(i)+nz-1
                    ! Ignore lower triangle entries
                    nz_ignored = 0
                    do j=this_lbound,this_ubound
                        if(this%ja(j)>=i) exit
                        nz_ignored = nz_ignored+1
                    enddo
                    ! Assign columns
                    A_II%irp(i+1)                                = A_XX_ubound-nz_ignored+1
                    A_II%ja (A_XX_lbound:A_XX_ubound-nz_ignored) = this%ja(this_lbound+nz_ignored:this_ubound)
                    A_II%nnz = A_II%nnz - nz_ignored + nz
                else
                    A_II%irp(i+1) = A_II%irp(i)
                endif

            else if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
                do j=this%irp(i), this%irp(i)+nz-1
                    ! Add upper triangle elements
                    k = this%ja(j)
                    A_II%ja(A_II%irp(i)) = k
                    A_II%irp(i) = A_II%irp(i)+1
                    if(i==k) cycle
                    ! Add lower triangle elements
                    A_II%ja(A_II%irp(k)) = i
                    A_II%irp(k) = A_II%irp(k)+1
                enddo
            endif

            ! Number of nnz of A_IG in row i
            nz = this%irp(i+1)-this%irp(i)-nz_offset 
            if(nz>0) then
                ! Calculate bounds
                A_XX_lbound = A_IG%irp(i);           A_XX_ubound = A_IG%irp(i)+nz-1
                this_lbound = this%irp(i)+nz_offset; this_ubound = this%irp(i+1)-1
                ! Assign columns
                A_IG%irp(i+1)                     = A_XX_ubound+1
                A_IG%ja (A_XX_lbound:A_XX_ubound) = this%ja (this_lbound:this_ubound)-num_col
                A_IG%nnz = A_IG%nnz + nz
            else
                A_IG%irp(i+1) = A_IG%irp(i)
            endif
        enddo

        if(this%get_symmetric_storage() .and. .not. A_II%get_symmetric_storage()) then
            do i=num_row, 1, -1
                A_II%irp(i+1) = A_II%irp(i)
            enddo
            A_II%irp(1) = 1
            A_II%nnz = A_II%irp(num_row+1)-1
        endif

        call memrealloc(A_II%nnz, A_II%ja, __FILE__, __LINE__)
        call memrealloc(A_IG%nnz, A_IG%ja, __FILE__, __LINE__)

        if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
            ! If input matrix has symmetric storage and one or more of the output diagonal matrix has not symmetric 
            A_GG%irp = 0
            ! Count total number of
            do i=num_row+1, total_rows
                do j=this%irp(i), this%irp(i+1)-1
                    k = this%ja(j)-num_col
                    if(k<0) cycle
                    A_GG%irp(i-num_row+1) = A_GG%irp(i-num_row+1)+1
                    if(i-num_row == k) cycle ! diagonal elements
                    A_GG%irp(k+1) = A_GG%irp(k+1)+1
                enddo
            enddo
            ! Build CSR header from counter
            A_GG%irp(1) = 1
            do i=1, A_GG%get_num_rows()
                A_GG%irp(i+1) = A_GG%irp(i+1)+A_GG%irp(i)
            enddo
        endif

        ! Loop (num_row:this%num_rows) to get A_GI and A_GG
        do i=num_row+1, total_rows
            nz_offset = 0
            ! Count the number of columns less than or equal to num_row
            do j=this%irp(i), this%irp(i+1)-1
                if(this%ja(j)> num_col) exit
                nz_offset = nz_offset + 1
            enddo

            ! Number of nnz of A_GI in row i-num_row
            nz = nz_offset 
            if(present(A_GI)) then
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_GI%irp(i-num_row); A_XX_ubound = A_GI%irp(i-num_row)+nz-1
                    this_lbound = this%irp(i);         this_ubound = this%irp(i)+nz-1
                    ! Assign columns
                    A_GI%irp(i-num_row+1)             = A_XX_ubound+1
                    A_GI%ja (A_XX_lbound:A_XX_ubound) = this%ja (this_lbound:this_ubound)
                    A_GI%nnz = A_GI%nnz + nz
                else
                    A_GI%irp(i-num_row+1) = A_GI%irp(i-num_row)
                endif
            endif

            ! Number of nnz of A_GG in row i-num_row
            nz = this%irp(i+1)-this%irp(i)-nz_offset 
            if(this%get_symmetric_storage() .eqv. A_GG%get_symmetric_storage()) then
                ! Number of nnz of A_GG in row i-num_row
                nz = this%irp(i+1)-this%irp(i)-nz_offset 
                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_GG%irp(i-num_row);   A_XX_ubound = A_GG%irp(i-num_row)+nz-1
                    this_lbound = this%irp(i)+nz_offset; this_ubound = this%irp(i+1)-1
                    ! Assign columns
                    A_GG%irp(i-num_row+1)             = A_XX_ubound+1
                    A_GG%ja (A_XX_lbound:A_XX_ubound) = this%ja (this_lbound:this_ubound)-num_col
                    A_GG%nnz = A_GG%nnz + nz
                else
                    A_GG%irp(i-num_row+1) = A_GG%irp(i-num_row)
                endif

            else if(.not. this%get_symmetric_storage() .and. A_GG%get_symmetric_storage()) then

                if(nz>0) then
                    ! Calculate bounds
                    A_XX_lbound = A_GG%irp(i-num_row);   A_XX_ubound = A_GG%irp(i-num_row)+nz-1
                    this_lbound = this%irp(i)+nz_offset; this_ubound = this%irp(i+1)-1
                    ! Ignore lower triangle entries
                    nz_ignored = 0
                    do j=this_lbound,this_ubound
                        if(this%ja(j)>=i) exit
                        nz_ignored = nz_ignored+1
                    enddo
                    ! Assign columns
                    A_GG%irp(i-num_row+1)             = A_XX_ubound-nz_ignored+1
                    A_GG%ja (A_XX_lbound:A_XX_ubound-nz_ignored) = this%ja (this_lbound+nz_ignored:this_ubound)-num_col
                    A_GG%nnz = A_GG%nnz - nz_ignored + nz
                else
                    A_GG%irp(i-num_row+1) = A_GG%irp(i-num_row)
                endif

            else if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
                do j=this%irp(i), this%irp(i)+nz-1
                    ! Add upper triangle elements
                    k = this%ja(j)-num_col
                    A_GG%ja(A_GG%irp(i-num_row)) = k
                    A_GG%irp(i-num_row) = A_GG%irp(i-num_row)+1
                    if(i-num_row==k) cycle
                    ! Add lower triangle elements
                    A_GG%ja(A_GG%irp(k)) = i-num_row
                    A_GG%irp(k) = A_GG%irp(k)+1
                enddo
            endif
        enddo

        if(this%get_symmetric_storage() .and. .not. A_GG%get_symmetric_storage()) then
            do i=A_GG%get_num_rows(), 1, -1
                A_GG%irp(i+1) = A_GG%irp(i)
            enddo
            A_GG%irp(1) = 1
            A_GG%nnz = A_GG%irp(A_GG%get_num_rows()+1)-1
        endif

        if(present(A_GI)) call memrealloc(A_GI%nnz, A_GI%ja, __FILE__, __LINE__)
        call memrealloc(A_GG%nnz, A_GG%ja, __FILE__, __LINE__)

        call A_II%set_state_assembled_symbolic()
        call A_IG%set_state_assembled_symbolic()
        if(present(A_GI)) call A_GI%set_state_assembled_symbolic()
        call A_GG%set_state_assembled_symbolic()
    end subroutine csr_sparse_matrix_split_2x2_symbolic_body


    subroutine csr_sparse_matrix_permute_and_split_2x2_numeric(this, num_row, num_col, perm, iperm, A_CC, A_CR, A_RC, A_RR)
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 and permute some columns and rows 
    !< given 2 permutation arrays (perm and iperm)
    !< 
    !< A = [A_CC A_RC]
    !<     [A_CR A_RR]
    !<
    !< this routine computes A_CC, A_RC, A_CR and A_RR given the global 
    !< matrix A. Note that A_CC, A_RC, A_CR are dense and A_RR is sparse
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        integer(ip),                           intent(in)    :: perm(:)
        integer(ip),                           intent(in)    :: iperm(:)
        real(rp),    allocatable,              intent(out)   :: A_CC(:,:)
        real(rp),    allocatable,              intent(out)   :: A_CR(:,:)
        real(rp),    allocatable,              intent(out)   :: A_RC(:,:)
        class(base_sparse_matrix_t),           intent(inout) :: A_RR
    !-----------------------------------------------------------------
        select type (A_RR)
            type is (csr_sparse_matrix_t)
                call this%permute_and_split_2x2_numeric_body(num_row, num_col, perm, iperm,  A_CC, A_CR, A_RC, A_RR)
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_permute_and_split_2x2_numeric


    subroutine csr_sparse_matrix_permute_and_split_2x2_numeric_body(this, num_row, num_col, perm, iperm, A_CC, A_CR, A_RC, A_RR) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 and permute some columns and rows 
    !< given 2 permutation arrays (perm and iperm)
    !< 
    !< A = [A_CC A_RC]
    !<     [A_CR A_RR]
    !<
    !< this routine computes A_CC, A_RC, A_CR and A_RR given the global 
    !< matrix A. Note that A_CC, A_RC, A_CR are dense and A_RR is sparse
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        integer(ip),                           intent(in)    :: perm(:)
        integer(ip),                           intent(in)    :: iperm(:)
        real(rp),    allocatable,              intent(out)   :: A_CC(:,:)
        real(rp),    allocatable,              intent(out)   :: A_CR(:,:)
        real(rp),    allocatable,              intent(out)   :: A_RC(:,:)
        class(csr_sparse_matrix_t),            intent(inout) :: A_RR
        logical                                              :: THIS_has_symmetric_storage
        logical                                              :: A_RR_has_symmetric_storage
        logical                                              :: symmetric
        integer(ip)                                          :: sign
        integer(ip)                                          :: total_num_rows
        integer(ip)                                          :: total_num_cols
        integer(ip)                                          :: A_RR_num_rows
        integer(ip)                                          :: A_RR_num_cols
        integer(ip)                                          :: perm_size
        integer(ip)                                          :: current_row
        integer(ip)                                          :: current_col
        integer(ip)                                          :: current_row_offset
        integer(ip)                                          :: permuted_row
        integer(ip)                                          :: permuted_col
        integer(ip)                                          :: permuted_row_offset
        integer(ip)                                          :: next_permuted_row_offset
        integer(ip)                                          :: nz_per_row
        integer(ip), allocatable                             :: link_list(:)
        integer(ip)                                          :: iret,i,j
    !-----------------------------------------------------------------
        assert(this%state_is_assembled() .or. this%state_is_assembled_symbolic())
        assert(A_RR%state_is_properties_set() .or. A_RR%state_is_assembled_symbolic())

        total_num_rows = this%get_num_rows()
        total_num_cols = this%get_num_cols()
        assert(num_row>0 .and. num_row<total_num_rows .and. num_col>0 .and. num_col<total_num_cols)
        perm_size = max(total_num_rows, total_num_cols)

        assert(size(perm) == perm_size)
        assert(size(iperm) == perm_size)

        ! Allocate A_CC, A_CR and A_RC arrays
        call memalloc(perm_size+2, link_list, __FILE__, __LINE__)
        call memalloc(num_row, num_col, A_CC, __FILE__, __LINE__);                A_CC = 0._rp
        call memalloc(num_row, total_num_cols-num_col, A_CR, __FILE__, __LINE__); A_CR = 0._rp
        call memalloc(total_num_rows-num_row, num_col, A_RC, __FILE__, __LINE__); A_RC = 0._rp

        A_RR_num_rows = total_num_rows-num_row
        A_RR_num_cols = total_num_cols-num_col
        THIS_has_symmetric_storage = this%get_symmetric_storage()
        A_RR_has_symmetric_storage = A_RR%get_symmetric_storage()
        if(A_RR%state_is_properties_set()) then
            call A_RR%set_num_rows(A_RR_num_rows)
            call A_RR%set_num_cols(A_RR_num_cols)
            if(THIS_has_symmetric_storage .and. .not. A_RR_has_symmetric_storage) then
                ! All diagonal elements in the original matrix must appear in the sparsity pattern
                call A_RR%allocate_numeric(2*this%get_nnz()-total_num_rows)
            else if(.not. THIS_has_symmetric_storage .and. A_RR_has_symmetric_storage) then
                ! All diagonal elements in the original matrix must appear in the sparsity pattern
                call A_RR%allocate_numeric((this%get_nnz()-total_num_rows)/2+total_num_rows)
            else
                call A_RR%allocate_numeric(this%get_nnz())
            endif
        else
            assert(A_RR%get_num_rows() == A_RR_num_rows .and. A_RR%get_num_cols() == A_RR_num_rows)
            call A_RR%allocate_values_body(A_RR%get_nnz())
        endif

        A_RR%irp = 0

        if(THIS_has_symmetric_storage .and. .not. A_RR_has_symmetric_storage) then
            ! Count nnz per row for A_RR matrix (transpose non diagonal elements)
            do permuted_row=num_row+1, total_num_rows
                current_row = iperm(permuted_row)
                do i = this%irp(current_row), this%irp(current_row+1)-1
                    permuted_col = perm(this%ja(i))
                    A_RR%Irp(permuted_row-num_row+1) = A_RR%irp(permuted_row-num_row+1)+1
                    if(permuted_row-num_row == permuted_col-num_col .or. permuted_col<=num_col) cycle
                    A_RR%irp(permuted_col-num_col+1) = A_RR%irp(permuted_col-num_col+1)+1
                enddo
            enddo

            ! Rebuild CSR header from nz counter
            A_RR%irp(1) = 1
            do permuted_row=1, total_num_rows-num_row
                A_RR%irp(permuted_row+1) = A_RR%irp(permuted_row+1)+A_RR%irp(permuted_row)
            enddo
        endif

        ! Loop on permuted rows to fill A_CC and A_CR arrays
        A_RR%irp(1) = 1
        do permuted_row=1, num_row
            current_row = iperm(permuted_row)
            current_row_offset = this%irp(current_row)
            do i = current_row_offset, this%irp(current_row+1)-1
                current_col = this%ja(i)
                permuted_col = perm(current_col)
                if ( permuted_col <= num_col ) then
                    A_CC(permuted_row,permuted_col) = this%val(i)
                    if(THIS_has_symmetric_storage .and. permuted_col<=num_row .and. &
                            permuted_row<=num_col) &
                            A_CC(permuted_col,permuted_row) = this%val(i)
                else
                    A_CR(permuted_row,permuted_col-num_col) = this%val(i)
                    if(THIS_has_symmetric_storage .and. permuted_col>num_row .and. &
                            permuted_row<=num_col) &
                            A_RC(permuted_col-num_col,permuted_row) = this%val(i)
                endif
            end do
        enddo

        ! Loop on permuted rows to fill A_RC array and A_RR sparse matrix
        do permuted_row=num_row+1, total_num_rows
            current_row = iperm(permuted_row)
            current_row_offset = this%irp(current_row)
            nz_per_row = 0

            if(THIS_has_symmetric_storage .and. .not. A_RR_has_symmetric_storage) then
                permuted_row_offset = A_RR%irp(permuted_row-num_row)
                ! Add permuted_cols to the permuted_row
                do i = current_row_offset, this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    permuted_col = perm(current_col)
                    if ( permuted_col <= num_col ) then
                        A_RC(permuted_row-num_row,permuted_col) = this%val(i)
                        if(THIS_has_symmetric_storage .and. (permuted_col<num_row) .and. &
                            (permuted_row>num_col)) &
                            A_CR(permuted_col,permuted_row-num_row) = this%val(i)
                    else
                        ! Add permuted_col to the permuted_row of A_RR
                        j = A_RR%irp(permuted_row-num_row)
                        A_RR%ja(j) = permuted_col-num_col
                        A_RR%val(j) = this%val(i)
                        A_RR%irp(permuted_row-num_row) = A_RR%irp(permuted_row-num_row)+1
                        if(permuted_row-num_row == permuted_col-num_col) cycle
                        ! Unsymmetrize element if is not in the diagonal
                        j = A_RR%irp(permuted_col-num_col)
                        A_RR%ja(j) = permuted_row-num_row
                        A_RR%val(j) = this%val(i)
                        A_RR%irp(permuted_col-num_col) = A_RR%irp(permuted_col-num_col)+1
                    endif
                enddo

            else
                permuted_row_offset = A_RR%irp(permuted_row-num_row)
                ! Add permuted_cols to the permuted_row
                do i = current_row_offset, this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    permuted_col = perm(current_col)
                    if ( permuted_col <= num_col ) then
                        A_RC(permuted_row-num_row,permuted_col) = this%val(i)
                        if(THIS_has_symmetric_storage .and. (permuted_col<num_row) .and. &
                            (permuted_row>num_col)) &
                            A_CR(permuted_col,permuted_row-num_row) = this%val(i)
                    else
                        if (A_RR_has_symmetric_storage .and. permuted_row-num_row>permuted_col-num_col) cycle
                        A_RR%ja(permuted_row_offset+nz_per_row) = permuted_col-num_col
                        A_RR%val(permuted_row_offset+nz_per_row) = this%val(i)
                        nz_per_row = nz_per_row + 1
                    endif
                end do
                A_RR%irp(permuted_row-num_row+1) = permuted_row_offset+nz_per_row

                ! Sort permuted columns of A_RR sparse matrix for each row
                if(nz_per_row>0) then
                    next_permuted_row_offset = permuted_row_offset+nz_per_row
                    call mergesort_link_list(nz_per_row,                         &
                        A_RR%ja(permuted_row_offset:next_permuted_row_offset-1), &
                        link_list,                                               &
                        iret)
                    if(iret == 0) call reorder_ip_rp_from_link_list(nz_per_row,   &
                        A_RR%val(permuted_row_offset:next_permuted_row_offset-1), &
                        A_RR%ja(permuted_row_offset:next_permuted_row_offset-1),  &
                        link_list)
                endif
            endif
        end do

        if(THIS_has_symmetric_storage .and. .not. A_RR_has_symmetric_storage) then
            ! In this case A_RR%irp is not well builded we have to revert it
            ! to get a valid CSR header
            do permuted_row=A_RR_num_rows, 1, -1
                A_RR%irp(permuted_row+1) = A_RR%irp(permuted_row)
            enddo
            A_RR%irp(1) = 1

            ! Sort permuted_cols for each permuted_row
            do permuted_row=1, A_RR_num_rows
                permuted_row_offset = A_RR%irp(permuted_row)
                next_permuted_row_offset = A_RR%irp(permuted_row+1)-1
                nz_per_row = next_permuted_row_offset-permuted_row_offset+1
                ! Sort permuted columns of A_RR sparse matrix for each row
                call mergesort_link_list(nz_per_row,                         &
                    A_RR%ja(permuted_row_offset:next_permuted_row_offset),   &
                    link_list,                                               &
                    iret)
                if(iret == 0) call reorder_ip_rp_from_link_list(nz_per_row,  &
                    A_RR%val(permuted_row_offset:next_permuted_row_offset),  &
                    A_RR%ja(permuted_row_offset:next_permuted_row_offset),   &
                    link_list)
            enddo
        endif

        A_RR%nnz = A_RR%irp(A_RR_num_rows+1)-1
        call memfree(link_list, __FILE__, __LINE__)

        ! Adjust A_RR columns and vals size to A_RR%nnz
        call memrealloc(A_RR%nnz, A_RR%ja, __FILE__, __LINE__)
        call memrealloc(A_RR%nnz, A_RR%val, __FILE__, __LINE__)

        call A_RR%set_state_assembled()
    end subroutine csr_sparse_matrix_permute_and_split_2x2_numeric_body


    subroutine csr_sparse_matrix_permute_and_split_2x2_symbolic(this, num_row, num_col, perm, iperm, A_RR) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 and permute some columns and rows
    !< given 2 permutation arrays (perm and iperm)
    !< 
    !< A = [A_CC A_RC]
    !<     [A_CR A_RR]
    !<
    !< this routine computes A_RR from the global matrix A
    !< A_CC, ACR and A_RC sparsity pattern calculation is not
    !< performed because they are dense matrices
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        integer(ip),                           intent(in)    :: perm(:)
        integer(ip),                           intent(in)    :: iperm(:)
        class(base_sparse_matrix_t),           intent(inout) :: A_RR
    !-----------------------------------------------------------------
        select type (A_RR)
            type is (csr_sparse_matrix_t)
                call this%permute_and_split_2x2_symbolic_body(num_row, num_col, perm, iperm, A_RR)
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_permute_and_split_2x2_symbolic


    subroutine csr_sparse_matrix_permute_and_split_2x2_symbolic_body(this, num_row, num_col, perm, iperm, A_RR) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 and permute some columns and rows
    !< given 2 permutation arrays (perm and iperm)
    !< 
    !< A = [A_CC A_RC]
    !<     [A_CR A_RR]
    !<
    !< this routine computes A_RR from the global matrix A
    !< A_CC, ACR and A_RC sparsity pattern calculation is not
    !< performed because they are dense matrices
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        integer(ip),                           intent(in)    :: perm(:)
        integer(ip),                           intent(in)    :: iperm(:)
        class(csr_sparse_matrix_t),            intent(inout) :: A_RR
        logical                                              :: THIS_has_symmetric_storage
        logical                                              :: A_RR_has_symmetric_storage
        integer(ip)                                          :: total_num_rows
        integer(ip)                                          :: total_num_cols
        integer(ip)                                          :: A_RR_num_rows
        integer(ip)                                          :: A_RR_num_cols
        integer(ip)                                          :: perm_size
        integer(ip)                                          :: current_row
        integer(ip)                                          :: current_col
        integer(ip)                                          :: current_row_offset
        integer(ip)                                          :: permuted_row
        integer(ip)                                          :: permuted_col
        integer(ip)                                          :: permuted_row_offset
        integer(ip)                                          :: next_permuted_row_offset
        integer(ip)                                          :: nz_per_row
        integer(ip), allocatable                             :: link_list(:)
        integer(ip)                                          :: iret,i,j
    !-----------------------------------------------------------------
        assert(this%state_is_assembled() .or. this%state_is_assembled_symbolic())
        assert(A_RR%state_is_properties_set())

        total_num_rows = this%get_num_rows()
        total_num_cols = this%get_num_cols()
        assert(num_row>0 .and. num_row<total_num_rows .and. num_col>0 .and. num_col<total_num_cols)
        perm_size = max(total_num_rows, total_num_cols)

        assert(size(perm) == perm_size)
        assert(size(iperm) == perm_size)

        A_RR_num_rows = total_num_rows-num_row
        A_RR_num_cols = total_num_cols-num_col
        call A_RR%set_num_rows(A_RR_num_rows)
        call A_RR%set_num_cols(A_RR_num_cols)
        THIS_has_symmetric_storage = this%get_symmetric_storage()
        A_RR_has_symmetric_storage = A_RR%get_symmetric_storage()

        call memalloc(perm_size+2, link_list, __FILE__, __LINE__)

        if(THIS_has_symmetric_storage .and. .not. A_RR_has_symmetric_storage) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            call A_RR%allocate_symbolic(2*this%get_nnz()-total_num_rows)
        else if(.not. THIS_has_symmetric_storage .and. A_RR_has_symmetric_storage) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            call A_RR%allocate_symbolic((this%get_nnz()-total_num_rows)/2+total_num_rows)
        else
            call A_RR%allocate_symbolic(this%get_nnz())
        endif

        A_RR%irp = 0

        if(THIS_has_symmetric_storage .and. .not. A_RR_has_symmetric_storage) then
            ! Count nnz per row for A_RR matrix (transpose non diagonal elements)
            do permuted_row=1, A_RR_num_rows
                current_row = iperm(permuted_row+num_row)
                do i = this%irp(current_row), this%irp(current_row+1)-1
                    permuted_col = perm(this%ja(i))-num_col
                    A_RR%Irp(permuted_row+1) = A_RR%irp(permuted_row+1)+1
                    if(permuted_row == permuted_col) cycle
                    A_RR%irp(permuted_col+1) = A_RR%irp(permuted_col+1)+1
                enddo
            enddo

            ! Build CSR header from counter
            A_RR%irp(1) = 1
            do permuted_row=1, A_RR_num_rows
                A_RR%irp(permuted_row+1) = A_RR%irp(permuted_row+1)+A_RR%irp(permuted_row)
            enddo

            do permuted_row=1, A_RR_num_rows
                current_row = iperm(permuted_row+num_row)
                current_row_offset = this%irp(current_row)
                nz_per_row = 0

                permuted_row_offset = A_RR%irp(permuted_row)
                ! Add permuted_cols to the permuted_row
                do i = current_row_offset, this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    permuted_col = perm(current_col)-num_col
                    if ( permuted_col > 0 ) then
                        ! Add permuted_col to the permuted_row of A_RR
                        j = A_RR%irp(permuted_row)
                        A_RR%ja(j) = permuted_col
                        A_RR%irp(permuted_row) = A_RR%irp(permuted_row)+1
                        if(permuted_row == permuted_col) cycle
                        ! Unsymmetrize element if is not in the diagonal
                        j = A_RR%irp(permuted_col)
                        A_RR%ja(j) = permuted_row
                        A_RR%irp(permuted_col) = A_RR%irp(permuted_col)+1
                    endif
                enddo
            enddo

            ! In this case A_RR%irp is not well builded we have to revert it
            ! to get a valid CSR header
            do permuted_row=A_RR_num_rows, 1, -1
                A_RR%irp(permuted_row+1) = A_RR%irp(permuted_row)
            enddo
            A_RR%irp(1) = 1

            ! Sort permuted_cols for each permuted_row
            do permuted_row=1, A_RR_num_rows
                permuted_row_offset = A_RR%irp(permuted_row)
                next_permuted_row_offset = A_RR%irp(permuted_row+1)-1
                nz_per_row = next_permuted_row_offset-permuted_row_offset+1
                ! Sort permuted columns of A_RR sparse matrix for each row
                call mergesort_link_list(nz_per_row,                         &
                    A_RR%ja(permuted_row_offset:next_permuted_row_offset), &
                    link_list,                                               &
                    iret)
                if(iret == 0) call reorder_ip_from_link_list(nz_per_row,     &
                    A_RR%ja(permuted_row_offset:next_permuted_row_offset), &
                    link_list)
            enddo

        else
            ! Loop on permuted rows to fill A_RR sparse matrix
            A_RR%irp(1) = 1
            do permuted_row=1, A_RR_num_rows
                current_row = iperm(permuted_row+num_row)
                current_row_offset = this%irp(current_row)
                permuted_row_offset = A_RR%irp(permuted_row)
                nz_per_row = 0

                ! Add permuted_cols to the permuted_row
                do i = current_row_offset, this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    permuted_col = perm(current_col) - num_col
                    if (permuted_col<1 .or. (A_RR_has_symmetric_storage .and. permuted_row>permuted_col)) cycle
                    A_RR%ja(permuted_row_offset+nz_per_row) = permuted_col
                    nz_per_row = nz_per_row + 1
                end do

                ! Sort permuted columns of A_RR sparse matrix for each row
                if(nz_per_row>0) then
                    next_permuted_row_offset = permuted_row_offset+nz_per_row
                    call mergesort_link_list(nz_per_row,                         &
                        A_RR%ja(permuted_row_offset:next_permuted_row_offset-1), &
                        link_list,                                               &
                        iret)
                    if(iret == 0) call reorder_ip_from_link_list(nz_per_row,   &
                        A_RR%ja(permuted_row_offset:next_permuted_row_offset-1),  &
                        link_list)
                endif
                A_RR%irp(permuted_row+1) = permuted_row_offset+nz_per_row
            end do

        endif

        A_RR%nnz = A_RR%irp(A_RR_num_rows+1)-1

        call memfree(link_list, __FILE__, __LINE__)

        ! Adjust A_RR columns size to A_RR%nnz
        call memrealloc(A_RR%nnz, A_RR%ja, __FILE__, __LINE__)

        call A_RR%set_state_assembled_symbolic()
    end subroutine csr_sparse_matrix_permute_and_split_2x2_symbolic_body


    subroutine csr_sparse_matrix_expand_matrix_numeric_array(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, C_T_val, I_nz, I_ia, I_ja, I_val, to)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),      intent(in)    :: this
        integer,                         intent(in)    :: C_T_num_cols
        integer,                         intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        real(rp),                        intent(in)    :: C_T_val(C_T_nz)
        integer,                         intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        real(rp),                        intent(in)    :: I_val(C_T_nz)
        class(base_sparse_matrix_t),     intent(inout) :: to
    !-----------------------------------------------------------------
        select type (to)
            type is (csr_sparse_matrix_t)
                call this%expand_matrix_numeric_array_body(C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, C_T_val, I_nz, I_ia, I_ja, I_val, to)
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_expand_matrix_numeric_array


    subroutine csr_sparse_matrix_expand_matrix_numeric_array_body(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, C_T_val, I_nz, I_ia, I_ja, I_val, to)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO format
    !< A = [A C_T]
    !<     [C  I ]
    !< Some considerations:
    !<  - C = transpose(C_T)
    !<  - I is a square matrix
    !<  - THIS (input) sparse matrix must be in ASSEMBLED state
    !<  - TO (output) sparse matrix must be in PROPERTIES_SET state
    !<  - C_T coordinate arrays (C_T_ia, C_T_ja and C_T_val) must 
    !<    have the same size (C_T_nz)
    !<  - I coordinate arrays (I_ia, I_ja and I_val) must 
    !<    have the same size (I_nz)
    !<  - Row index arrays (X_ia) must be in ascendent order
    !<  - Column index arrays (X_ja) must be in ascendent order for 
    !<    each row
    !<  - For each C_T row index (C_T_ia): 1<=C_T_ia(i)<=this%get_num_rows()
    !<  - For each C_T column index (C_T_ja): 1<=C_T_ja(i)<=C_T_num_cols
    !<  - For each I row and column index (I_ia and I_ja): 
    !<    1<=I_ia(i) and I_ia(i)<=C_T_num_cols
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),      intent(in)    :: this
        integer,                         intent(in)    :: C_T_num_cols
        integer,                         intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        real(rp),                        intent(in)    :: C_T_val(C_T_nz)
        integer,                         intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        real(rp),                        intent(in)    :: I_val(I_nz)
        class(csr_sparse_matrix_t),      intent(inout) :: to
        integer                                        :: f, i, j, k
        integer                                        :: initial_num_rows
        integer                                        :: initial_num_cols
        integer                                        :: C_T_num_rows
        integer                                        :: previous_ia
        integer                                        :: previous_ja
        integer                                        :: new_nz
        integer                                        :: nz_per_row
        integer                                        :: current_nz_per_row
        integer                                        :: nz_offset
        integer                                        :: last_visited_row
        integer                                        :: current_row
        integer                                        :: current_col
        integer                                        :: next_row
        integer                                        :: next_row_offset
        integer                                        :: last_visited_row_offset
        integer                                        :: C_irp(C_T_num_cols)
        integer                                        :: I_irp(C_T_num_cols)
        integer                                        :: I_counter
        integer                                        :: nz_per_row_counter(C_T_num_cols)
        integer                                        :: nz_counter
        logical                                        :: sorted
        logical                                        :: symmetric_storage
    !-----------------------------------------------------------------
        assert(this%state_is_assembled())
        assert(to%state_is_properties_set() .or. to%state_is_assembled_symbolic())

        initial_num_rows = this%get_num_rows()
        initial_num_cols = this%get_num_cols()

    !-----------------------------------------------------------------
    ! Set properties to the expanded matrix
    !-----------------------------------------------------------------

        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        symmetric_storage = to%get_symmetric_storage()

#ifdef DEBUG
        C_T_num_rows     = maxval(C_T_ia)
        assert(C_T_num_rows == initial_num_rows .and. initial_num_rows>0 .and. initial_num_cols>0)

    !-----------------------------------------------------------------
    ! Check if (C_T) ia and ja arrays are sorted by rows
    ! It also counts number or colums per row for C matrix
    !-----------------------------------------------------------------
        sorted = .true.
        previous_ia = 1
        previous_ja = 0
        do i=1, C_T_nz
            if(previous_ia /= C_T_ia(i)) previous_ja = 0
            if((C_T_ia(i)>initial_num_rows) .or. (C_T_ja(i)>C_T_num_cols) .or. &
               (previous_ia>C_T_ia(i)) .or. (previous_ja>=C_T_ja(i))) then
                sorted = .false.
                exit
            endif
            previous_ia = C_T_ia(i)
            previous_ja = C_T_ja(i)
        enddo
        assert(sorted)
    !-----------------------------------------------------------------
    ! Check if (I) ia and ja arrays are sorted by rows
    ! It also counts number or colums per row for I matrix
    !-----------------------------------------------------------------
        previous_ia = 1
        previous_ja = 0
        do i=1, I_nz
            if(previous_ia /= I_ia(i)) previous_ja = 0
            if((I_ia(i)>C_T_num_cols) .or. (I_ja(i)>C_T_num_cols) .or. &
               (previous_ia>I_ia(i)) .or. (previous_ja>=I_ja(i))) then
                sorted = .false.
                exit
            endif
            previous_ia = I_ia(i)
            previous_ja = I_ja(i)
        enddo
        assert(sorted)
#endif

        ! Count number of elements of I finally added to the output matrix
        I_irp = 0
        do i=1, I_nz
            if(symmetric_storage .and. I_ia(i)>I_ja(i)) cycle
            I_irp(I_ia(i)) = I_irp(I_ia(i)) + 1
        enddo

    !-----------------------------------------------------------------
    ! Alloc to%irp with the new number of rows and to%ja with the new number of nnz
    !-----------------------------------------------------------------
        if(           this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=max(this%nnz, 2*this%nnz-initial_num_rows)+2*C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=(this%nnz-initial_num_rows)/2+initial_num_rows+C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            new_nz=this%nnz+2*C_T_nz+sum(I_irp)
        else if(      this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            new_nz=this%nnz+C_T_nz+sum(I_irp)
        endif
        if(to%state_is_properties_set()) then
            call memalloc(initial_num_rows+C_T_num_cols+1, to%irp, __FILE__, __LINE__)
            call memalloc(new_nz, to%ja, __FILE__, __LINE__)
        else
            assert(size(to%irp) == initial_num_rows+C_T_num_cols+1)
            assert(size(to%ja) == new_nz)
        endif
        call memalloc(new_nz, to%val, __FILE__, __LINE__)

    !-----------------------------------------------------------------
    ! Expand  C_T matrix (Add columns to existing rows)
    !-----------------------------------------------------------------
        nz_counter = 0
        C_Irp = 0
        if(this%get_symmetric_storage() .eqv. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(:initial_num_rows) = this%irp(:initial_num_rows)

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            nz_counter = 0
            current_row = C_T_ia(1)
            if(current_row/=1) then
                nz_counter = this%irp(current_row)-this%irp(1)
                to%ja(1:nz_counter) = this%ja(this%irp(1):this%irp(current_row)-1)
                to%val(1:nz_counter) = this%val(this%irp(1):this%irp(current_row)-1)
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            current_nz_per_row = 0
            do i=1,C_T_nz
                nz_per_row = nz_per_row+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 
                current_nz_per_row      = next_row_offset-last_visited_row_offset             ! Number of nnz per row before expanding the matrix
                ! Append existing columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+current_nz_per_row) = this%ja(last_visited_row_offset:next_row_offset-1)
                to%val(nz_counter+1:nz_counter+current_nz_per_row) = this%val(last_visited_row_offset:next_row_offset-1)
                nz_counter = nz_counter+current_nz_per_row
                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                to%val(nz_counter+1:nz_counter+nz_per_row) = C_T_val(i-nz_per_row+1:i)
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row:initial_num_rows) = to%irp(next_row:initial_num_rows) + nz_per_row
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            to%nnz = this%nnz + C_T_nz
            ! If the last visited row is not the last row append the rest of the ja values
            if(last_visited_row/=initial_num_rows+1) then
                to%ja(nz_counter+1:to%nnz) = this%ja(this%irp(last_visited_row):this%irp(initial_num_rows+1)-1)
                to%val(nz_counter+1:to%nnz) = this%val(this%irp(last_visited_row):this%irp(initial_num_rows+1)-1)
            endif
        else if (.not. this%get_symmetric_storage() .and. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(1) = 1

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            ! Filter lower triangle entries
            current_row = C_T_ia(1)
            if(current_row/=1) then
                do i=1, current_row-1
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    to%val(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%val(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            do i=1,C_T_nz
                nz_per_row = nz_per_row+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 

                ! Append existing columns into the current row of the expanded matrix
                ! Filter lower triangle entries
                do j=last_visited_row, current_row
                    do nz_offset=this%irp(j), this%irp(j+1)-1
                        if(this%ja(nz_offset)>=j) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(j+1)-nz_offset) = this%ja(nz_offset:this%irp(j+1)-1)
                    to%val(nz_counter+1:nz_counter+this%irp(j+1)-nz_offset) = this%val(nz_offset:this%irp(j+1)-1)
                    nz_counter = nz_counter + this%irp(j+1)-nz_offset
                    to%irp(j+1) = nz_counter+1
                enddo

                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                to%val(nz_counter+1:nz_counter+nz_per_row) = C_T_val(i-nz_per_row+1:i)
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row) = nz_counter+1
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            ! If the last visited row is not the last row append the rest of the ja values
            ! Filter lower triangle entries
            if(last_visited_row/=initial_num_rows+1) then
                do i=last_visited_row, initial_num_rows
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    to%val(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%val(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            to%nnz = nz_counter
        else if (this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            to%irp = 0; to%irp(1) = 1
    !-----------------------------------------------------------------        
    ! Count number of elements in each row for the 
    ! lower (new) + upper triangle, till initial_num_rows
    !-----------------------------------------------------------------        
            do current_row=1, initial_num_rows
                ! lower triangle, transposed elements (new)
                do i=this%irp(current_row), this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    if(current_row==current_col) cycle
                    to%irp(current_col+1) = to%irp(current_col+1)+1
                enddo
                ! Upper triangle elements
                to%irp(current_row+1) = to%irp(current_row+1)+this%irp(current_row+1)-this%irp(current_row)
            enddo

            ! Add the number of elements Of C_T in each row
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                to%irp(current_row+1) = to%irp(current_row+1)+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
            enddo

            ! Convert counter to CSR header
            do current_row=1, max(initial_num_rows, initial_num_cols)
                to%irp(current_row+1) = to%irp(current_row+1)+to%irp(current_row)
            enddo

    !-----------------------------------------------------------------        
    ! Add elements to JA and VAL arrays
    !-----------------------------------------------------------------        
            ! Add elements from the original matrix
            do current_row=1, initial_num_rows
                do i=this%irp(current_row),this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    to%ja(to%irp(current_row)) = current_col
                    to%val(to%irp(current_row)) = this%val(i)
                    to%irp(current_row) = to%irp(current_row)+1
                    if(current_row == current_col) cycle
                    to%ja(to%irp(current_col)) = current_row
                    to%val(to%irp(current_col)) = this%val(i)
                    to%irp(current_col) = to%irp(current_col)+1
                enddo
            enddo

            ! Add C_T_nz elements
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                current_col = C_T_ja(i)
                to%ja(to%irp(current_row)) = initial_num_cols+ current_col
                to%val(to%irp(current_row)) = C_T_val(i)
                to%irp(current_row) = to%irp(current_row)+1
            enddo

            ! Recalculate CSR irp header
            to%nnz = to%irp(initial_num_rows)-1
            do current_row=initial_num_rows, 1, -1
                to%irp(current_row+1)=to%irp(current_row)
            enddo
            to%irp(1) = 1
        endif
        to%irp(initial_num_rows+1) = to%nnz+1

    !-----------------------------------------------------------------
    ! Loop to expand with C  and I matrices (Append new rows)
    !-----------------------------------------------------------------        
        nz_per_row_counter = 0
        I_counter = 1
        do i=1,C_T_num_cols
            current_row = initial_num_rows+i+1
            if(symmetric_storage) then
                ! If symmetric_storage, only upper triangle of I matrix will be appended
                nz_per_row = I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        if(I_ja(I_counter)>=I_ia(I_counter)) then
                            to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                            to%val(to%irp(current_row-1)+nz_per_row_counter(i)) = I_val(I_counter)
                            nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        endif
                        I_counter = I_counter+1
                        if(I_counter>I_nz) exit
                    enddo
                endif
            else
                ! If not symmetric_storage, both, C and I matrix will be appended
                nz_per_row = C_irp(i)
                ! C_T_ja are the rows of C
                ! C_T_ia are the cols of C
                do j=1,C_T_nz
                    if(C_T_ja(j)==i) then
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = C_T_ia(j)
                        to%val(to%irp(current_row-1)+nz_per_row_counter(i)) = C_T_val(j)
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        if(nz_per_row_counter(i)>=nz_per_row) exit
                    endif
                enddo
                nz_per_row = nz_per_row + I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                        to%val(to%irp(current_row-1)+nz_per_row_counter(i)) = I_val(I_counter)
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        I_counter = I_counter+1
                        if(I_counter>I_nz .or. nz_per_row_counter(i)>=nz_per_row) exit
                    enddo
                endif
            endif
            to%irp(current_row) = to%irp(current_row-1)+nz_per_row_counter(i)
            nz_per_row = 0
        enddo
    !-----------------------------------------------------------------    
    ! Update matrix properties
    !-----------------------------------------------------------------    
        to%nnz = to%nnz + sum(nz_per_row_counter)
        to%irp(initial_num_rows+C_T_num_cols+1) = to%nnz+1
        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        call to%set_state_assembled()
    end subroutine csr_sparse_matrix_expand_matrix_numeric_array_body


    subroutine csr_sparse_matrix_expand_matrix_numeric_coo(this, C_T, to, I)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),          intent(in)    :: this
        type(coo_sparse_matrix_t),           intent(in)    :: C_T
        class(base_sparse_matrix_t),         intent(inout) :: to
        type(coo_sparse_matrix_t), optional, intent(in)    :: I
    !-----------------------------------------------------------------
        select type (to)
            type is (csr_sparse_matrix_t)
                call this%expand_matrix_numeric_coo_body(C_T, to, I)
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_expand_matrix_numeric_coo


    subroutine csr_sparse_matrix_expand_matrix_numeric_coo_body(this, C_T_coo, to, I_coo)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO format
    !< A = [A C_T]
    !<     [C  I ]
    !< Some considerations:
    !<  - C = transpose(C_T)
    !<  - I is a square matrix
    !<  - THIS (input) sparse matrix must be in ASSEMBLED state
    !<  - TO (output) sparse matrix must be in PROPERTIES_SET state
    !<  - C_T is a COO sparse matrix assembled and sorted by rows
    !<  - I is a square COO sparse matrix assembled and sorted by rows
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),        intent(in)  :: this
        type(coo_sparse_matrix_t), target, intent(in)  :: C_T_coo
        class(csr_sparse_matrix_t),      intent(inout) :: to
        type(coo_sparse_matrix_t), optional, target, intent(in)  :: I_coo
        integer(ip)                                    :: C_T_num_cols
        integer(ip)                                    :: C_T_num_rows
        integer(ip)                                    :: C_T_nz
        integer(ip), pointer                           :: C_T_ia(:)
        integer(ip), pointer                           :: C_T_ja(:)
        real(rp),    pointer                           :: C_T_val(:)
        integer(ip)                                    :: I_nz
        integer(ip), pointer                           :: I_ia(:)
        integer(ip), pointer                           :: I_ja(:)
        real(rp),    pointer                           :: I_val(:)
        integer(ip)                                    :: f, i, j, k
        integer(ip)                                    :: initial_num_rows
        integer(ip)                                    :: initial_num_cols
        integer(ip)                                    :: previous_ia
        integer(ip)                                    :: previous_ja
        integer(ip)                                    :: new_nz
        integer(ip)                                    :: nz_per_row
        integer(ip)                                    :: current_nz_per_row
        integer(ip)                                    :: nz_offset
        integer(ip)                                    :: last_visited_row
        integer(ip)                                    :: current_row
        integer(ip)                                    :: current_col
        integer(ip)                                    :: next_row
        integer(ip)                                    :: next_row_offset
        integer(ip)                                    :: last_visited_row_offset
        integer(ip), allocatable                       :: C_irp(:)
        integer(ip), allocatable                       :: I_irp(:)
        integer(ip)                                    :: I_counter
        integer(ip), allocatable                       :: nz_per_row_counter(:)
        integer(ip)                                    :: nz_counter
        logical(ip)                                    :: sorted
        logical(ip)                                    :: symmetric_storage
    !-----------------------------------------------------------------
        assert(this%state_is_assembled())
        assert(C_T_coo%state_is_assembled() .and. C_T_coo%is_by_rows())
        if(present(I_coo)) then
            assert(I_coo%state_is_assembled() .and. I_coo%is_by_rows())
            assert(C_T_coo%get_num_cols() == I_coo%get_num_cols() .and. I_coo%get_num_rows() == I_coo%get_num_cols())
        endif
        assert(to%state_is_properties_set() .or. to%state_is_assembled_symbolic())

        C_T_num_cols =  C_T_coo%get_num_cols()
        C_T_num_rows =  C_T_coo%get_num_rows()
        C_T_nz       =  C_T_coo%get_nnz()
        C_T_ia       => C_T_coo%ia
        C_T_ja       => C_T_coo%ja
        C_T_val      => C_T_coo%val
        if(present(I_coo)) then
            I_nz         =  I_coo%get_nnz()
            I_ia         => I_coo%ia
            I_ja         => I_coo%ja
            I_val        => I_coo%val
        else
            I_nz = C_T_num_cols
            allocate(I_ia(C_T_num_cols+1))
            allocate(I_ja(C_T_num_cols))
            allocate(I_val(C_T_num_cols))
            I_val = 0.0_rp
            do i=1, I_nz
                I_ia(i) = i
                I_ja(i) = i
            enddo
            I_ia(I_nz+1) = I_nz+1
        endif

        initial_num_rows = this%get_num_rows()
        initial_num_cols = this%get_num_cols()
        assert(C_T_num_rows == initial_num_rows .and. initial_num_rows>0 .and. initial_num_cols>0)

        allocate(C_irp(C_T_num_cols))
        allocate(I_irp(C_T_num_cols))
        allocate(nz_per_row_counter(C_T_num_cols))
    !-----------------------------------------------------------------
    ! Set properties to the expanded matrix
    !-----------------------------------------------------------------

        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        symmetric_storage = to%get_symmetric_storage()

        ! Count number of elements of I finally added to the output matrix
        I_irp = 0
        do i=1, I_nz
            if(symmetric_storage .and. I_ia(i)>I_ja(i)) cycle
            I_irp(I_ia(i)) = I_irp(I_ia(i)) + 1
        enddo

        
    !-----------------------------------------------------------------
    ! Alloc to%irp with the new number of rows and to%ja with the new number of nnz
    !-----------------------------------------------------------------
            
        if(           this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=max(this%nnz, 2*this%nnz-initial_num_rows)+2*C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=(this%nnz-initial_num_rows)/2+initial_num_rows+C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            new_nz=this%nnz+2*C_T_nz+sum(I_irp)
        else if(      this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            new_nz=this%nnz+C_T_nz+sum(I_irp)
        endif
        if(to%state_is_properties_set()) then
            call memalloc(initial_num_rows+C_T_num_cols+1, to%irp, __FILE__, __LINE__)
            call memalloc(new_nz, to%ja, __FILE__, __LINE__)
        else
            assert(size(to%irp) == initial_num_rows+C_T_num_cols+1)
            assert(size(to%ja) == new_nz)
        endif
        call memalloc(new_nz, to%val, __FILE__, __LINE__)

    !-----------------------------------------------------------------
    ! Expand  C_T matrix (Add columns to existing rows)
    !-----------------------------------------------------------------
        nz_counter = 0
        C_Irp = 0
        if(this%get_symmetric_storage() .eqv. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(:initial_num_rows) = this%irp(:initial_num_rows)

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            nz_counter = 0
            current_row = C_T_ia(1)
            if(current_row/=1) then
                nz_counter = this%irp(current_row)-this%irp(1)
                to%ja(1:nz_counter) = this%ja(this%irp(1):this%irp(current_row)-1)
                to%val(1:nz_counter) = this%val(this%irp(1):this%irp(current_row)-1)
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            current_nz_per_row = 0
            do i=1,C_T_nz
                nz_per_row = nz_per_row+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 
                current_nz_per_row      = next_row_offset-last_visited_row_offset             ! Number of nnz per row before expanding the matrix
                ! Append existing columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+current_nz_per_row) = this%ja(last_visited_row_offset:next_row_offset-1)
                to%val(nz_counter+1:nz_counter+current_nz_per_row) = this%val(last_visited_row_offset:next_row_offset-1)
                nz_counter = nz_counter+current_nz_per_row
                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                to%val(nz_counter+1:nz_counter+nz_per_row) = C_T_val(i-nz_per_row+1:i)
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row:initial_num_rows) = to%irp(next_row:initial_num_rows) + nz_per_row
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            to%nnz = this%nnz + C_T_nz
            ! If the last visited row is not the last row append the rest of the ja values
            if(last_visited_row/=initial_num_rows+1) then
                to%ja(nz_counter+1:to%nnz) = this%ja(this%irp(last_visited_row):this%irp(initial_num_rows+1)-1)
                to%val(nz_counter+1:to%nnz) = this%val(this%irp(last_visited_row):this%irp(initial_num_rows+1)-1)
            endif
        else if (.not. this%get_symmetric_storage() .and. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(1) = 1

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            ! Filter lower triangle entries
            current_row = C_T_ia(1)
            if(current_row/=1) then
                do i=1, current_row-1
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    to%val(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%val(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            do i=1,C_T_nz
                nz_per_row = nz_per_row+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 

                ! Append existing columns into the current row of the expanded matrix
                ! Filter lower triangle entries
                do j=last_visited_row, current_row
                    do nz_offset=this%irp(j), this%irp(j+1)-1
                        if(this%ja(nz_offset)>=j) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(j+1)-nz_offset) = this%ja(nz_offset:this%irp(j+1)-1)
                    to%val(nz_counter+1:nz_counter+this%irp(j+1)-nz_offset) = this%val(nz_offset:this%irp(j+1)-1)
                    nz_counter = nz_counter + this%irp(j+1)-nz_offset
                    to%irp(j+1) = nz_counter+1
                enddo

                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                to%val(nz_counter+1:nz_counter+nz_per_row) = C_T_val(i-nz_per_row+1:i)
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row) = nz_counter+1
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            ! If the last visited row is not the last row append the rest of the ja values
            ! Filter lower triangle entries
            if(last_visited_row/=initial_num_rows+1) then
                do i=last_visited_row, initial_num_rows
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    to%val(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%val(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            to%nnz = nz_counter
        else if (this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            to%irp = 0; to%irp(1) = 1
    !-----------------------------------------------------------------        
    ! Count number of elements in each row for the 
    ! lower (new) + upper triangle, till initial_num_rows
    !-----------------------------------------------------------------        
            do current_row=1, initial_num_rows
                ! lower triangle, transposed elements (new)
                do i=this%irp(current_row), this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    if(current_row==current_col) cycle
                    to%irp(current_col+1) = to%irp(current_col+1)+1
                enddo
                ! Upper triangle elements
                to%irp(current_row+1) = to%irp(current_row+1)+this%irp(current_row+1)-this%irp(current_row)
            enddo

            ! Add the number of elements Of C_T in each row
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                to%irp(current_row+1) = to%irp(current_row+1)+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
            enddo

            ! Convert counter to CSR header
            do current_row=1, max(initial_num_rows, initial_num_cols)
                to%irp(current_row+1) = to%irp(current_row+1)+to%irp(current_row)
            enddo

    !-----------------------------------------------------------------        
    ! Add elements to JA and VAL arrays
    !-----------------------------------------------------------------        
            ! Add elements from the original matrix
            do current_row=1, initial_num_rows
                do i=this%irp(current_row),this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    to%ja(to%irp(current_row)) = current_col
                    to%val(to%irp(current_row)) = this%val(i)
                    to%irp(current_row) = to%irp(current_row)+1
                    if(current_row == current_col) cycle
                    to%ja(to%irp(current_col)) = current_row
                    to%val(to%irp(current_col)) = this%val(i)
                    to%irp(current_col) = to%irp(current_col)+1
                enddo
            enddo

            ! Add C_T_nz elements
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                current_col = C_T_ja(i)
                to%ja(to%irp(current_row)) = initial_num_cols+ current_col
                to%val(to%irp(current_row)) = C_T_val(i)
                to%irp(current_row) = to%irp(current_row)+1
            enddo

            ! Recalculate CSR irp header
            to%nnz = to%irp(initial_num_rows)-1
            do current_row=initial_num_rows, 1, -1
                to%irp(current_row+1)=to%irp(current_row)
            enddo
            to%irp(1) = 1
        endif
        to%irp(initial_num_rows+1) = to%nnz+1

    !-----------------------------------------------------------------
    ! Loop to expand with C  and I matrices (Append new rows)
    !-----------------------------------------------------------------        
        nz_per_row_counter = 0
        I_counter = 1
        do i=1,C_T_num_cols
            current_row = initial_num_rows+i+1
            if(symmetric_storage) then
                ! If symmetric_storage, only upper triangle of I matrix will be appended
                nz_per_row = I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        if(I_ja(I_counter)>=I_ia(I_counter)) then
                            to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                            to%val(to%irp(current_row-1)+nz_per_row_counter(i)) = I_val(I_counter)
                            nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        endif
                        I_counter = I_counter+1
                        if(I_counter>I_nz) exit
                    enddo
                endif
            else
                ! If not symmetric_storage, both, C and I matrix will be appended
                nz_per_row = C_irp(i)
                ! C_T_ja are the rows of C
                ! C_T_ia are the cols of C
                do j=1,C_T_nz
                    if(C_T_ja(j)==i) then
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = C_T_ia(j)
                        to%val(to%irp(current_row-1)+nz_per_row_counter(i)) = C_T_val(j)
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        if(nz_per_row_counter(i)>=nz_per_row) exit
                    endif
                enddo
                nz_per_row = nz_per_row + I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                        to%val(to%irp(current_row-1)+nz_per_row_counter(i)) = I_val(I_counter)
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        I_counter = I_counter+1
                        if(I_counter>I_nz .or. nz_per_row_counter(i)>=nz_per_row) exit
                    enddo
                endif
            endif
            to%irp(current_row) = to%irp(current_row-1)+nz_per_row_counter(i)
            nz_per_row = 0
        enddo
    !-----------------------------------------------------------------    
    ! Update matrix properties
    !-----------------------------------------------------------------    
        to%nnz = to%nnz + sum(nz_per_row_counter)
        to%irp(initial_num_rows+C_T_num_cols+1) = to%nnz+1
        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        call to%set_state_assembled()

        if(allocated(C_irp)) deallocate(C_irp)
        if(allocated(I_irp)) deallocate(I_irp)
        if(allocated(nz_per_row_counter)) deallocate(nz_per_row_counter)
        if(.not. present(I_coo)) then
            deallocate(I_ia)
            deallocate(I_ja)
            deallocate(I_val)
        endif
    end subroutine csr_sparse_matrix_expand_matrix_numeric_coo_body


    subroutine csr_sparse_matrix_expand_matrix_symbolic_array(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, I_nz, I_ia, I_ja, to)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),      intent(in)    :: this
        integer(ip),                     intent(in)    :: C_T_num_cols
        integer(ip),                     intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        integer(ip),                     intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        class(base_sparse_matrix_t),     intent(inout) :: to
    !-----------------------------------------------------------------
        select type (to)
            type is (csr_sparse_matrix_t)
                call this%expand_matrix_symbolic_array_body(C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, I_nz, I_ia, I_ja, to)
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_expand_matrix_symbolic_array


    subroutine csr_sparse_matrix_expand_matrix_symbolic_array_body(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, I_nz, I_ia, I_ja, to)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !< Some considerations:
    !<  - C = transpose(C_T)
    !<  - I is a square matrix
    !<  - THIS (input sparse matrix must be in ASSEMBLED or 
    !<    ASSEMBLED_SYMBOLIC state
    !<  - TO (output) sparse matrix must be in PROPERTIES_SET state
    !<  - C_T coordinate arrays (C_T_ia, C_T_ja and C_T_val) must 
    !<    have the same size (C_T_nz)
    !<  - I coordinate arrays (I_ia, I_ja and I_val) must 
    !<    have the same size (I_nz)
    !<  - Row index arrays (X_ia) must be in ascendent order
    !<  - Column index arrays (X_ja) must be in ascendent order for 
    !<    each row
    !<  - For each C_T row index (C_T_ia): 1<=C_T_ia(i)<=this%get_num_rows()
    !<  - For each C_T column index (C_T_ja): 1<=C_T_ja(i)<=C_T_num_cols
    !<  - For each I row and column index (I_ia and I_ja): 
    !<    1<=I_ia(i) and I_ia(i)<=C_T_num_cols
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),      intent(in)    :: this
        integer,                         intent(in)    :: C_T_num_cols
        integer,                         intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        integer,                         intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        class(csr_sparse_matrix_t),      intent(inout) :: to
        integer                                        :: f, i, j, k
        integer                                        :: initial_num_rows
        integer                                        :: initial_num_cols
        integer                                        :: C_T_num_rows
        integer                                        :: previous_ia
        integer                                        :: previous_ja
        integer                                        :: new_nz
        integer                                        :: nz_per_row
        integer                                        :: current_nz_per_row
        integer                                        :: nz_offset
        integer                                        :: last_visited_row
        integer                                        :: current_row
        integer                                        :: current_col
        integer                                        :: next_row
        integer                                        :: next_row_offset
        integer                                        :: last_visited_row_offset
        integer                                        :: C_irp(C_T_num_cols)
        integer                                        :: I_irp(C_T_num_cols)
        integer                                        :: I_counter
        integer                                        :: nz_per_row_counter(C_T_num_cols)
        integer                                        :: nz_counter
        logical                                        :: sorted
        logical                                        :: symmetric_storage
    !-----------------------------------------------------------------
        assert(this%state_is_assembled() .or. this%state_is_assembled_symbolic())
        assert(to%state_is_properties_set())

        initial_num_rows = this%get_num_rows()
        initial_num_cols = this%get_num_cols()

    !-----------------------------------------------------------------
    ! Set properties to the expanded matrix
    !-----------------------------------------------------------------
        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        symmetric_storage = to%get_symmetric_storage()
#ifdef DEBUG
        C_T_num_rows     = maxval(C_T_ia)
        assert(C_T_num_rows == initial_num_rows .and. initial_num_rows>0 .and. initial_num_cols>0)

    !-----------------------------------------------------------------
    ! Check if (C_T) ia and ja arrays are sorted by rows
    ! It also counts number or colums per row for C matrix
    !-----------------------------------------------------------------
        sorted = .true.
        previous_ia = 1
        previous_ja = 0
        do i=1, C_T_nz
            if(previous_ia /= C_T_ia(i)) previous_ja = 0
            if((C_T_ia(i)>initial_num_rows) .or. (C_T_ja(i)>C_T_num_cols) .or. &
               (previous_ia>C_T_ia(i)) .or. (previous_ja>=C_T_ja(i))) then
                sorted = .false.
                exit
            endif
            previous_ia = C_T_ia(i)
            previous_ja = C_T_ja(i)
        enddo
        assert(sorted)
    !-----------------------------------------------------------------
    ! Check if (I) ia and ja arrays are sorted by rows
    ! It also counts number or colums per row for I matrix
    !-----------------------------------------------------------------
        previous_ia = 1
        previous_ja = 0
        do i=1, I_nz
            if(previous_ia /= I_ia(i)) previous_ja = 0
            if((I_ia(i)>C_T_num_cols) .or. (I_ja(i)>C_T_num_cols) .or. &
               (previous_ia>I_ia(i)) .or. (previous_ja>=I_ja(i))) then
                sorted = .false.
                exit
            endif
            previous_ia = I_ia(i)
            previous_ja = I_ja(i)
        enddo
        assert(sorted)
#endif

        ! Count number of elements of I finally added to the output matrix
        I_irp = 0
        do i=1, I_nz
            if(symmetric_storage .and. I_ia(i)>I_ja(i)) cycle
            I_irp(I_ia(i)) = I_irp(I_ia(i)) + 1
        enddo

    !-----------------------------------------------------------------
    ! Alloc to%irp with the new number of rows and to%ja with the new number of nnz
    !-----------------------------------------------------------------
        call memalloc(initial_num_rows+C_T_num_cols+1, to%irp, __FILE__, __LINE__)
        if(           this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=max(this%nnz,2*this%nnz-initial_num_rows)+2*C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=(this%nnz-initial_num_rows)/2+initial_num_rows+C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            new_nz=this%nnz+2*C_T_nz+sum(I_irp)
        else if(      this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            new_nz=this%nnz+C_T_nz+sum(I_irp)
        endif
        call memalloc(new_nz, to%ja, __FILE__, __LINE__)

    !-----------------------------------------------------------------
    ! Expand  C_T matrix (Add columns to existing rows)
    !-----------------------------------------------------------------
        nz_counter = 0
        C_irp = 0
        if(this%get_symmetric_storage() .eqv. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(:initial_num_rows) = this%irp(:initial_num_rows)

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            current_row = C_T_ia(1)
            if(current_row/=1) then
                nz_counter = this%irp(current_row)-this%irp(1)
                to%ja(1:nz_counter) = this%ja(this%irp(1):this%irp(current_row)-1)
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            current_nz_per_row = 0
            do i=1,C_T_nz
                nz_per_row = nz_per_row+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 
                current_nz_per_row      = next_row_offset-last_visited_row_offset             ! Number of nnz per row before expanding the matrix
                ! Append existing columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+current_nz_per_row) = this%ja(last_visited_row_offset:next_row_offset-1)
                nz_counter = nz_counter+current_nz_per_row
                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row:initial_num_rows) = to%irp(next_row:initial_num_rows) + nz_per_row
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            ! If the last visited row is not the last row append the rest of the ja values
            to%nnz = this%nnz + C_T_nz
            if(last_visited_row/=initial_num_rows+1) then
                to%ja(nz_counter+1:to%nnz) = this%ja(this%irp(last_visited_row):this%irp(initial_num_rows+1)-1)
            endif
        else if (.not. this%get_symmetric_storage() .and. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(1) = 1

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            ! Filter lower triangle entries
            current_row = C_T_ia(1)
            if(current_row/=1) then
                do i=1, current_row-1
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            do i=1,C_T_nz
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                nz_per_row = nz_per_row+1
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 

                ! Append existing columns into the current row of the expanded matrix
                ! Filter lower triangle entries
                do j=last_visited_row, current_row
                    do nz_offset=this%irp(j), this%irp(j+1)-1
                        if(this%ja(nz_offset)>=j) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(j+1)-nz_offset) = this%ja(nz_offset:this%irp(j+1)-1)
                    nz_counter = nz_counter + this%irp(j+1)-nz_offset
                    to%irp(j+1) = nz_counter+1
                enddo

                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row) = nz_counter+1
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            ! If the last visited row is not the last row append the rest of the ja values
            ! Filter lower triangle entries
            if(last_visited_row/=initial_num_rows+1) then
                do i=last_visited_row, initial_num_rows
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            to%nnz = nz_counter
        else if (this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp = 0; to%irp(1) = 1
    !-----------------------------------------------------------------        
    ! Count number of elements in each row for the 
    ! lower (new) + upper triangle, till initial_num_rows
    !-----------------------------------------------------------------        
            do current_row=1, initial_num_rows
                ! lower triangle, transposed elements (new)
                do i=this%irp(current_row), this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    if(current_row==current_col) cycle
                    to%irp(current_col+1) = to%irp(current_col+1)+1
                enddo
                ! Upper triangle elements
                to%irp(current_row+1) = to%irp(current_row+1)+this%irp(current_row+1)-this%irp(current_row)
            enddo

            ! Add the number of elements Of C_T in each row
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                to%irp(current_row+1) = to%irp(current_row+1)+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
            enddo

            ! Convert counter to CSR header
            do current_row=1, max(initial_num_rows, initial_num_cols)
                to%irp(current_row+1) = to%irp(current_row+1)+to%irp(current_row)
            enddo

    !-----------------------------------------------------------------        
    ! Add elements to JA array
    !-----------------------------------------------------------------        
            ! Add elements from the original matrix
            do current_row=1, initial_num_rows
                do i=this%irp(current_row),this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    to%ja(to%irp(current_row)) = current_col
                    to%irp(current_row) = to%irp(current_row)+1
                    if(current_row == current_col) cycle
                    to%ja(to%irp(current_col)) = current_row
                    to%irp(current_col) = to%irp(current_col)+1
                enddo
            enddo

            ! Add C_T_nz elements
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                current_col = C_T_ja(i)
                to%ja(to%irp(current_row)) = initial_num_cols+ current_col
                to%irp(current_row) = to%irp(current_row)+1
            enddo

            ! Recalculate CSR irp header
            to%nnz = to%irp(initial_num_rows)-1
            do current_row=initial_num_rows, 1, -1
                to%irp(current_row+1)=to%irp(current_row)
            enddo
            to%irp(1) = 1

        endif
        to%irp(initial_num_rows+1) = to%nnz+1

    !-----------------------------------------------------------------
    ! Loop to expand with C  and I matrices (Append new rows)
    !-----------------------------------------------------------------
        nz_per_row_counter = 0
        I_counter = 1
        do i=1,C_T_num_cols
            current_row = initial_num_rows+i+1
            if(symmetric_storage) then
                ! If symmetric_storage, only upper triangle of I matrix will be appended
                nz_per_row = I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        if(I_ja(I_counter)>=I_ia(I_counter)) then
                            to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                            nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        endif
                        I_counter = I_counter+1
                        if(I_counter>I_nz) exit
                    enddo
                endif
            else
                ! If not symmetric_storage, both, C and I matrix will be appended
                nz_per_row = C_irp(i)
                ! C_T_ja are the rows of C
                ! C_T_ia are the cols of C
                do j=1,C_T_nz
                    if(C_T_ja(j)==i) then
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = C_T_ia(j)
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        if(nz_per_row_counter(i)>=nz_per_row) exit
                    endif
                enddo
                nz_per_row = nz_per_row + I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        I_counter = I_counter+1
                        if(I_counter>I_nz .or. nz_per_row_counter(i)>=nz_per_row) exit
                    enddo
                endif
            endif
            to%irp(current_row) = to%irp(current_row-1)+nz_per_row_counter(i)
            nz_per_row = 0
        enddo
    !-----------------------------------------------------------------
    ! Update matrix properties
    !-----------------------------------------------------------------
        to%nnz = to%nnz + sum(nz_per_row_counter)
        to%irp(initial_num_rows+C_T_num_cols+1) = to%nnz+1
        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        call to%set_state_assembled_symbolic()
    end subroutine csr_sparse_matrix_expand_matrix_symbolic_array_body


    subroutine csr_sparse_matrix_expand_matrix_symbolic_coo(this, C_T, to, I)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),          intent(in)    :: this
        type(coo_sparse_matrix_t),           intent(in)    :: C_T
        type(coo_sparse_matrix_t), optional, intent(in)    :: I
        class(base_sparse_matrix_t),         intent(inout) :: to
    !-----------------------------------------------------------------
        select type (to)
            type is (csr_sparse_matrix_t)
                call this%expand_matrix_symbolic_coo_body(C_T, to, I)
            class DEFAULT
                check(.false.)
        end select
    end subroutine csr_sparse_matrix_expand_matrix_symbolic_coo


    subroutine csr_sparse_matrix_expand_matrix_symbolic_coo_body(this, C_T_coo, to, I_coo)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO format
    !< A = [A C_T]
    !<     [C  I ]
    !< Some considerations:
    !<  - C = transpose(C_T)
    !<  - I is a square matrix
    !<  - THIS (input) sparse matrix must be in ASSEMBLED state
    !<  - TO (output) sparse matrix must be in PROPERTIES_SET state
    !<  - C_T is a COO sparse matrix assembled and sorted by rows
    !<  - I is a square COO sparse matrix assembled and sorted by rows
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),      intent(in)    :: this
        type(coo_sparse_matrix_t), target, intent(in)  :: C_T_coo
        class(csr_sparse_matrix_t),      intent(inout) :: to
        type(coo_sparse_matrix_t), optional, target, intent(in)  :: I_coo
        integer(ip)                                    :: C_T_num_cols
        integer(ip)                                    :: C_T_num_rows
        integer(ip)                                    :: C_T_nz
        integer(ip), pointer                           :: C_T_ia(:)
        integer(ip), pointer                           :: C_T_ja(:)
        integer(ip)                                    :: I_nz
        integer(ip), pointer                           :: I_ia(:)
        integer(ip), pointer                           :: I_ja(:)
        integer(ip)                                    :: f, i, j, k
        integer(ip)                                    :: initial_num_rows
        integer(ip)                                    :: initial_num_cols
        integer(ip)                                    :: previous_ia
        integer(ip)                                    :: previous_ja
        integer(ip)                                    :: new_nz
        integer(ip)                                    :: nz_per_row
        integer(ip)                                    :: current_nz_per_row
        integer(ip)                                    :: nz_offset
        integer(ip)                                    :: last_visited_row
        integer(ip)                                    :: current_row
        integer(ip)                                    :: current_col
        integer(ip)                                    :: next_row
        integer(ip)                                    :: next_row_offset
        integer(ip)                                    :: last_visited_row_offset
        integer(ip), allocatable                       :: C_irp(:)
        integer(ip), allocatable                       :: I_irp(:)
        integer(ip)                                    :: I_counter
        integer(ip), allocatable                       :: nz_per_row_counter(:)
        integer(ip)                                    :: nz_counter
        logical(ip)                                    :: sorted
        logical(ip)                                    :: symmetric_storage
    !-----------------------------------------------------------------
        assert(this%state_is_assembled() .or. this%state_is_assembled_symbolic())
        assert(to%state_is_properties_set())
        assert(C_T_coo%state_is_assembled() .or. C_T_coo%state_is_assembled_symbolic())
        assert(C_T_coo%is_by_rows())
        if(present(I_coo)) then
            assert(I_coo%state_is_assembled() .or. I_coo%state_is_assembled_symbolic())
            assert(I_coo%is_by_rows())
            assert(C_T_coo%get_num_cols() == I_coo%get_num_cols() .and. I_coo%get_num_rows() == I_coo%get_num_cols())
        endif

        C_T_num_cols =  C_T_coo%get_num_cols()
        C_T_num_rows =  C_T_coo%get_num_rows()
        C_T_nz       =  C_T_coo%get_nnz()
        C_T_ia       => C_T_coo%ia
        C_T_ja       => C_T_coo%ja
        if(present(I_coo)) then
            I_nz         =  I_coo%get_nnz()
            I_ia         => I_coo%ia
            I_ja         => I_coo%ja
        else
            I_nz = C_T_num_cols
            allocate(I_ia(C_T_num_cols+1))
            allocate(I_ja(C_T_num_cols))
            do i=1, I_nz
                I_ia(i) = i
                I_ja(i) = i
            enddo
            I_ia(I_nz+1) = I_nz+1
        endif

        initial_num_rows = this%get_num_rows()
        initial_num_cols = this%get_num_cols()
        assert(C_T_num_rows == initial_num_rows .and. initial_num_rows>0 .and. initial_num_cols>0)

        allocate(C_irp(C_T_num_cols))
        allocate(I_irp(C_T_num_cols))
        allocate(nz_per_row_counter(C_T_num_cols))
    !-----------------------------------------------------------------
    ! Set properties to the expanded matrix
    !-----------------------------------------------------------------
        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        symmetric_storage = to%get_symmetric_storage()

        ! Count number of elements of I finally added to the output matrix
        I_irp = 0
        do i=1, I_nz
            if(symmetric_storage .and. I_ia(i)>I_ja(i)) cycle
            I_irp(I_ia(i)) = I_irp(I_ia(i)) + 1
        enddo

    !-----------------------------------------------------------------
    ! Alloc to%irp with the new number of rows and to%ja with the new number of nnz
    !-----------------------------------------------------------------
        call memalloc(initial_num_rows+C_T_num_cols+1, to%irp, __FILE__, __LINE__)
        if(           this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=max(this%nnz,2*this%nnz-initial_num_rows)+2*C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            ! All diagonal elements in the original matrix must appear in the sparsity pattern
            new_nz=(this%nnz-initial_num_rows)/2+initial_num_rows+C_T_nz+sum(I_irp)
        else if(.not. this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            new_nz=this%nnz+2*C_T_nz+sum(I_irp)
        else if(      this%get_symmetric_storage() .and.       to%get_symmetric_storage()) then
            new_nz=this%nnz+C_T_nz+sum(I_irp)
        endif
        call memalloc(new_nz, to%ja, __FILE__, __LINE__)

    !-----------------------------------------------------------------
    ! Expand  C_T matrix (Add columns to existing rows)
    !-----------------------------------------------------------------
        nz_counter = 0
        C_irp = 0
        if(this%get_symmetric_storage() .eqv. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(:initial_num_rows) = this%irp(:initial_num_rows)

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            current_row = C_T_ia(1)
            if(current_row/=1) then
                nz_counter = this%irp(current_row)-this%irp(1)
                to%ja(1:nz_counter) = this%ja(this%irp(1):this%irp(current_row)-1)
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            current_nz_per_row = 0
            do i=1,C_T_nz
                nz_per_row = nz_per_row+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 
                current_nz_per_row      = next_row_offset-last_visited_row_offset             ! Number of nnz per row before expanding the matrix
                ! Append existing columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+current_nz_per_row) = this%ja(last_visited_row_offset:next_row_offset-1)
                nz_counter = nz_counter+current_nz_per_row
                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row:initial_num_rows) = to%irp(next_row:initial_num_rows) + nz_per_row
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            ! If the last visited row is not the last row append the rest of the ja values
            to%nnz = this%nnz + C_T_nz
            if(last_visited_row/=initial_num_rows+1) then
                to%ja(nz_counter+1:to%nnz) = this%ja(this%irp(last_visited_row):this%irp(initial_num_rows+1)-1)
            endif
        else if (.not. this%get_symmetric_storage() .and. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp(1) = 1

            ! If the current_row of C_T is different to 1: Copy the original ja from 1:current_row_offset
            ! Filter lower triangle entries
            current_row = C_T_ia(1)
            if(current_row/=1) then
                do i=1, current_row-1
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            last_visited_row = current_row

            ! Loop over C_T to expand the matrix
            nz_per_row = 0
            do i=1,C_T_nz
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
                nz_per_row = nz_per_row+1
                if(i/=C_T_nz) then
                    if(C_T_ia(i) == C_T_ia(i+1)) cycle ! count new zeros in the same row
                endif

                current_row             = C_T_ia(i)                                           ! current row or C_T
                next_row                = current_row+1                                       ! next row of C_T
                last_visited_row_offset = this%irp(last_visited_row)                          ! ja offset for the last visited row of C_T
                next_row_offset         = this%irp(next_row)                                  ! ja offset for the next row of C_T 

                ! Append existing columns into the current row of the expanded matrix
                ! Filter lower triangle entries
                do j=last_visited_row, current_row
                    do nz_offset=this%irp(j), this%irp(j+1)-1
                        if(this%ja(nz_offset)>=j) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(j+1)-nz_offset) = this%ja(nz_offset:this%irp(j+1)-1)
                    nz_counter = nz_counter + this%irp(j+1)-nz_offset
                    to%irp(j+1) = nz_counter+1
                enddo

                ! Append new columns into the current row of the expanded matrix
                to%ja(nz_counter+1:nz_counter+nz_per_row) = C_T_ja(i-nz_per_row+1:i) + initial_num_cols
                nz_counter = nz_counter + nz_per_row
                ! Add the new nnz to irp
                to%irp(next_row) = nz_counter+1
                nz_per_row = 0
                last_visited_row  = next_row
            enddo

            ! If the last visited row is not the last row append the rest of the ja values
            ! Filter lower triangle entries
            if(last_visited_row/=initial_num_rows+1) then
                do i=last_visited_row, initial_num_rows
                    do nz_offset=this%irp(i), this%irp(i+1)-1
                        if(this%ja(nz_offset)>=i) exit
                    enddo
                    to%ja(nz_counter+1:nz_counter+this%irp(i+1)-nz_offset) = this%ja(nz_offset:this%irp(i+1)-1)
                    nz_counter = nz_counter + this%irp(i+1)-nz_offset
                    to%irp(i+1) = nz_counter+1
                enddo
            endif
            to%nnz = nz_counter
        else if (this%get_symmetric_storage() .and. .not. to%get_symmetric_storage()) then
            ! Initialize irp
            to%irp = 0; to%irp(1) = 1
    !-----------------------------------------------------------------        
    ! Count number of elements in each row for the 
    ! lower (new) + upper triangle, till initial_num_rows
    !-----------------------------------------------------------------        
            do current_row=1, initial_num_rows
                ! lower triangle, transposed elements (new)
                do i=this%irp(current_row), this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    if(current_row==current_col) cycle
                    to%irp(current_col+1) = to%irp(current_col+1)+1
                enddo
                ! Upper triangle elements
                to%irp(current_row+1) = to%irp(current_row+1)+this%irp(current_row+1)-this%irp(current_row)
            enddo

            ! Add the number of elements Of C_T in each row
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                to%irp(current_row+1) = to%irp(current_row+1)+1
                C_irp(C_T_ja(i)) = C_irp(C_T_ja(i)) + 1 ! Count C rows (transposed  C_T)
            enddo

            ! Convert counter to CSR header
            do current_row=1, max(initial_num_rows, initial_num_cols)
                to%irp(current_row+1) = to%irp(current_row+1)+to%irp(current_row)
            enddo

    !-----------------------------------------------------------------        
    ! Add elements to JA array
    !-----------------------------------------------------------------        
            ! Add elements from the original matrix
            do current_row=1, initial_num_rows
                do i=this%irp(current_row),this%irp(current_row+1)-1
                    current_col = this%ja(i)
                    to%ja(to%irp(current_row)) = current_col
                    to%irp(current_row) = to%irp(current_row)+1
                    if(current_row == current_col) cycle
                    to%ja(to%irp(current_col)) = current_row
                    to%irp(current_col) = to%irp(current_col)+1
                enddo
            enddo

            ! Add C_T_nz elements
            do i=1, C_T_nz
                current_row = C_T_ia(i)
                current_col = C_T_ja(i)
                to%ja(to%irp(current_row)) = initial_num_cols+ current_col
                to%irp(current_row) = to%irp(current_row)+1
            enddo

            ! Recalculate CSR irp header
            to%nnz = to%irp(initial_num_rows)-1
            do current_row=initial_num_rows, 1, -1
                to%irp(current_row+1)=to%irp(current_row)
            enddo
            to%irp(1) = 1

        endif
        to%irp(initial_num_rows+1) = to%nnz+1

    !-----------------------------------------------------------------
    ! Loop to expand with C  and I matrices (Append new rows)
    !-----------------------------------------------------------------
        nz_per_row_counter = 0
        I_counter = 1
        do i=1,C_T_num_cols
            current_row = initial_num_rows+i+1
            if(symmetric_storage) then
                ! If symmetric_storage, only upper triangle of I matrix will be appended
                nz_per_row = I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        if(I_ja(I_counter)>=I_ia(I_counter)) then
                            to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                            nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        endif
                        I_counter = I_counter+1
                        if(I_counter>I_nz) exit
                    enddo
                endif
            else
                ! If not symmetric_storage, both, C and I matrix will be appended
                nz_per_row = C_irp(i)
                ! C_T_ja are the rows of C
                ! C_T_ia are the cols of C
                do j=1,C_T_nz
                    if(C_T_ja(j)==i) then
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = C_T_ia(j)
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        if(nz_per_row_counter(i)>=nz_per_row) exit
                    endif
                enddo
                nz_per_row = nz_per_row + I_irp(i)

                if(I_counter<=I_nz) then
                    do while (I_ia(I_counter)==i)
                        to%ja(to%irp(current_row-1)+nz_per_row_counter(i)) = I_ja(I_counter)+initial_num_cols
                        nz_per_row_counter(i) = nz_per_row_counter(i) + 1
                        I_counter = I_counter+1
                        if(I_counter>I_nz .or. nz_per_row_counter(i)>=nz_per_row) exit
                    enddo
                endif
            endif
            to%irp(current_row) = to%irp(current_row-1)+nz_per_row_counter(i)
            nz_per_row = 0
        enddo
    !-----------------------------------------------------------------
    ! Update matrix properties
    !-----------------------------------------------------------------
        to%nnz = to%nnz + sum(nz_per_row_counter)
        to%irp(initial_num_rows+C_T_num_cols+1) = to%nnz+1
        call to%set_num_rows(initial_num_rows+C_T_num_cols)
        call to%set_num_cols(initial_num_cols+C_T_num_cols)
        call to%set_state_assembled_symbolic()

        if(allocated(C_irp)) deallocate(C_irp)
        if(allocated(I_irp)) deallocate(I_irp)
        if(allocated(nz_per_row_counter)) deallocate(nz_per_row_counter)
        if(.not. present(I_coo)) then
            deallocate(I_ia)
            deallocate(I_ja)
        endif
    end subroutine csr_sparse_matrix_expand_matrix_symbolic_coo_body


    subroutine csr_sparse_matrix_extract_diagonal(this, diagonal)
    !-----------------------------------------------------------------
    !< Return the diagonal of a CSR sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in)    :: this
        real(rp),   allocatable,    intent(inout) :: diagonal(:)
        integer(ip)                               :: diagonal_size
        integer(ip)                               :: row
        integer(ip)                               :: row_start_offset
        integer(ip)                               :: row_end_offset
        integer(ip)                               :: col_offset_in_row
    !-----------------------------------------------------------------
        assert(this%state_is_assembled())
        if(allocated(diagonal)) deallocate(diagonal)
        diagonal_size = min(this%get_num_rows(), this%get_num_cols())
        allocate(diagonal(diagonal_size))

        do row=1, diagonal_size
            row_start_offset = this%irp(row)
            row_end_offset = this%irp(row+1)-1
            col_offset_in_row = binary_search(row,row_end_offset-row_start_offset+1,this%ja(row_start_offset:row_end_offset))
            if (col_offset_in_row==-1) then
                diagonal(row) = 0.0_rp
            else
                diagonal(row) = this%val(row_start_offset+col_offset_in_row-1)
            endif
        enddo
    end subroutine csr_sparse_matrix_extract_diagonal


    subroutine csr_sparse_matrix_free_coords(this)
    !-----------------------------------------------------------------
    !< Clean coords of CSR sparse matrix format derived type
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout)  :: this
    !-----------------------------------------------------------------
        if(allocated(this%irp)) call memfree (this%irp, __FILE__, __LINE__)
        if(allocated(this%ja))  call memfree (this%ja,  __FILE__, __LINE__)
    end subroutine csr_sparse_matrix_free_coords


    subroutine csr_sparse_matrix_free_val(this)
    !-----------------------------------------------------------------
    !< Free values of CSR sparse matrix format derived type
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout)  :: this
    !-----------------------------------------------------------------
        if(allocated(this%val)) call memfree (this%val, __FILE__, __LINE__)
    end subroutine csr_sparse_matrix_free_val


    subroutine csr_sparse_matrix_print(this,lunou, only_graph)
    !-----------------------------------------------------------------
    !< Print a CSR matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),  intent(in) :: this
        integer(ip),                 intent(in) :: lunou
        logical, optional,           intent(in) :: only_graph
        logical                                 :: print_vals
        integer(ip)                             :: i,j
    !-----------------------------------------------------------------
        print_vals = .true.; if(present(only_graph)) print_vals = .not. only_graph
        write (lunou, '(a)')     '********************************************'
        write (lunou, '(a)')     '************* CSR data structure ***********'
        write (lunou, '(a)')     '********************************************'
        write (lunou, '(a,i10)') 'Number of rows:', this%get_num_rows()
        write (lunou, '(a,i10)') 'Number of cols:', this%get_num_cols()
        write (lunou, '(a,i10)') 'Number of non zeros (nnz):', this%get_nnz()
        write (lunou, '(a,i2)')  'Sign:', this%get_sign()
        write (lunou, '(a,l2)')  'Symmetric:', this%is_symmetric()
        write (lunou, '(a,l2)')  'Symmetric storage:', this%get_symmetric_storage()
    
    
        write (lunou, '(a)')     'Rows list (irp):'
        if(allocated(this%irp)) then
            write (lunou, *)    this%irp(1:this%get_num_rows()+1)
        else
            write (lunou,'(A)') 'Not allocated'
        endif
    
        write (lunou, '(a)')      'Columns list (ja):'
        if(allocated(this%ja)) then
            write (lunou, *)    this%ja(1:this%get_nnz())
        else
            write (lunou,'(A)') 'Not allocated'
        endif

        if(print_vals) then
            write (lunou, '(a)')      'Values list (val):'
            if(allocated(this%val)) then
                write (lunou, *)    this%val(1:this%get_nnz())
            else
                write (lunou,'(A)') 'Not allocated'
            endif
        endif
    end subroutine csr_sparse_matrix_print


    subroutine csr_sparse_matrix_print_matrix_market_body (this, lunou, ng, l2g)
        class(csr_sparse_matrix_t), intent(in) :: this
        integer(ip),                intent(in) :: lunou
        integer(ip), optional,      intent(in) :: ng
        integer(ip), optional,      intent(in) :: l2g (*)
        integer(ip) :: i, j
        integer(ip) :: nr, nc
    
        if ( present(ng) ) then 
            nr = ng
            nc = ng
        else
            nr = this%get_num_rows()
            nc = this%get_num_cols()
        end if

        write (lunou,'(a)') '%%MatrixMarket matrix coordinate real general'
            write (lunou,*) nr,nc,this%irp(this%get_num_rows()+1)-1
            do i=1,this%get_num_rows()
                do j=this%irp(i),this%irp(i+1)-1
                    if (present(l2g)) then
                        write(lunou,'(i12, i12, e32.25)') l2g(i), l2g(this%ja(j)), this%val(j)
                    else
                        write(lunou,'(i12, i12, e32.25)') i, this%ja(j), this%val(j)
                    end if
                end do
            end do


    end subroutine csr_sparse_matrix_print_matrix_market_body

    subroutine csr_sparse_matrix_create_iterator(this,iterator)
      class(csr_sparse_matrix_t)          , target     , intent(in)    :: this
      class(base_sparse_matrix_iterator_t), allocatable, intent(inout) :: iterator
      
      type(csr_sparse_matrix_iterator_t), allocatable :: csr_iterator

      allocate(csr_iterator)
      call csr_iterator%create(this)
      call move_alloc(from=csr_iterator, to=iterator)
    end subroutine csr_sparse_matrix_create_iterator

    function csr_sparse_matrix_get_entry(this, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Get the value in the entry (ia,ja) in the sparse matrix
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in)  :: this
        integer(ip),                intent(in)  :: ia
        integer(ip),                intent(in)  :: ja
        real(rp),                   intent(out) :: val
        logical                                 :: csr_sparse_matrix_get_entry

        integer(ip)                             :: ipaux,i1,i2,nr,nc
    !-----------------------------------------------------------------
        ! Ignore out of bounds entries
        if (ia<1 .or. ia>this%get_num_rows() .or. ja<1 .or.  ja>this%get_num_cols() .or. &
             (this%get_symmetric_storage() .and. ia>ja)) then
           csr_sparse_matrix_get_entry = .false.
        end if

        i1 = this%irp(ia)
        i2 = this%irp(ia+1)
        nc = i2-i1
        ipaux = binary_search(ja,nc,this%ja(i1:i2-1))

        if (ipaux>0) then
           val=this%val(i1+ipaux-1)
           csr_sparse_matrix_get_entry = .true.
        else
           csr_sparse_matrix_get_entry = .false.
        end if
      end function csr_sparse_matrix_get_entry
      
    !---------------------------------------------------------------------
    !< CSR_SPARSE_MATRIX_ITERATOR PROCEDURES
    !---------------------------------------------------------------------
    subroutine csr_sparse_matrix_iterator_create(this,csr_matrix)
      class(csr_sparse_matrix_iterator_t), intent(inout) :: this
      type(csr_sparse_matrix_t) , target, intent(in)    :: csr_matrix
      this%matrix => csr_matrix
      call this%init()
    end subroutine csr_sparse_matrix_iterator_create

    subroutine csr_sparse_matrix_iterator_init(this)
      class(csr_sparse_matrix_iterator_t), intent(inout) :: this
      
      this%row_index = 1
      this%nnz_index = 1
    end subroutine csr_sparse_matrix_iterator_init

    subroutine csr_sparse_matrix_iterator_free(this)
      class(csr_sparse_matrix_iterator_t), intent(inout) :: this
      
      this%row_index = -1
      this%nnz_index = -1
      this%matrix => NULL()
    end subroutine csr_sparse_matrix_iterator_free

    subroutine csr_sparse_matrix_iterator_next(this)
      class(csr_sparse_matrix_iterator_t), intent(inout) :: this

      assert( this%nnz_index < this%matrix%irp(this%row_index+1))

      this%nnz_index = this%nnz_index + 1
      if (this%nnz_index == this%matrix%irp(this%row_index+1)) then
         this%row_index = this%row_index + 1
      end if     
    end subroutine csr_sparse_matrix_iterator_next

    function csr_sparse_matrix_iterator_has_finished(this)
      class(csr_sparse_matrix_iterator_t), intent(in) :: this
      logical :: csr_sparse_matrix_iterator_has_finished
      
      csr_sparse_matrix_iterator_has_finished = (this%nnz_index > this%matrix%nnz)
    end function csr_sparse_matrix_iterator_has_finished

    function csr_sparse_matrix_iterator_get_row(this)
      class(csr_sparse_matrix_iterator_t), intent(in) :: this
      integer(ip) :: csr_sparse_matrix_iterator_get_row

      csr_sparse_matrix_iterator_get_row = this%row_index
    end function csr_sparse_matrix_iterator_get_row

    function csr_sparse_matrix_iterator_get_column(this)
      class(csr_sparse_matrix_iterator_t), intent(in) :: this
      integer(ip) :: csr_sparse_matrix_iterator_get_column

      csr_sparse_matrix_iterator_get_column = this%matrix%ja(this%nnz_index)
    end function csr_sparse_matrix_iterator_get_column

    function csr_sparse_matrix_iterator_get_entry(this)
      class(csr_sparse_matrix_iterator_t), intent(in) :: this
      real(rp) :: csr_sparse_matrix_iterator_get_entry

      csr_sparse_matrix_iterator_get_entry = this%matrix%val(this%nnz_index)
    end function csr_sparse_matrix_iterator_get_entry

    subroutine csr_sparse_matrix_iterator_set_entry(this,new_value)
      class(csr_sparse_matrix_iterator_t), intent(in) :: this
      real(rp)                           , intent(in) :: new_value

      this%matrix%val(this%nnz_index) = new_value
    end subroutine csr_sparse_matrix_iterator_set_entry

end module csr_sparse_matrix_names

