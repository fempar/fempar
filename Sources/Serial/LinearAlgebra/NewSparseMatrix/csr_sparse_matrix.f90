module csr_sparse_matrix_names

USE types_names
USE memor_names
USE vector_names
USE serial_scalar_array_names
USE base_sparse_matrix_names

implicit none

# include "debug.i90"

private


    type, extends(base_sparse_matrix_t) :: csr_sparse_matrix_t
        integer(ip), private       :: nnz = 0                     !< Number of non zeros
        integer(ip), allocatable   :: irp(:)                      !< Row pointers
        integer(ip), allocatable   :: ja(:)                       !< Column indices        
        real(rp),    allocatable   :: val(:)                      !< Values
    contains
    private
        procedure, public :: is_by_rows                              => csr_sparse_matrix_is_by_rows
        procedure, public :: is_by_cols                              => csr_sparse_matrix_is_by_cols
        procedure, public :: set_nnz                                 => csr_sparse_matrix_set_nnz
        procedure, public :: get_nnz                                 => csr_sparse_matrix_get_nnz
        procedure, public :: copy_to_coo                             => csr_sparse_matrix_copy_to_coo
        procedure, public :: copy_from_coo                           => csr_sparse_matrix_copy_from_coo
        procedure, public :: move_to_coo                             => csr_sparse_matrix_move_to_coo
        procedure, public :: move_from_coo                           => csr_sparse_matrix_move_from_coo
        procedure, public :: move_to_fmt                             => csr_sparse_matrix_move_to_fmt
        procedure, public :: move_from_fmt                           => csr_sparse_matrix_move_from_fmt
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
        procedure, public :: split_2x2                               => csr_sparse_matrix_split_2x2
        procedure, public :: free_coords                             => csr_sparse_matrix_free_coords
        procedure, public :: free_val                                => csr_sparse_matrix_free_val
        procedure, public :: apply_body                              => csr_sparse_matrix_apply_body
        procedure, public :: print_matrix_market_body                => csr_sparse_matrix_print_matrix_market_body
        procedure, public :: print                                   => csr_sparse_matrix_print
    end type csr_sparse_matrix_t

public :: csr_sparse_matrix_t

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


    subroutine csr_sparse_matrix_copy_to_coo(this, to)
    !-----------------------------------------------------------------
    !< Copy this (CSR) -> to (COO)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in)    :: this
        class(coo_sparse_matrix_t), intent(inout) :: to
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
    end subroutine csr_sparse_matrix_copy_to_coo


    subroutine csr_sparse_matrix_copy_from_coo(this, from)
    !-----------------------------------------------------------------
    !< Copy from (COO) -> this (CSR)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        class(coo_sparse_matrix_t), intent(in)    :: from
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
            return      
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
        call this%set_state(from%get_state())
        call memfree(itmp, __FILE__, __LINE__)
    end subroutine csr_sparse_matrix_copy_from_coo


    subroutine csr_sparse_matrix_move_to_coo(this, to)
    !-----------------------------------------------------------------
    !< Move this (CSR) -> to (COO)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        class(coo_sparse_matrix_t), intent(inout) :: to
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
    end subroutine csr_sparse_matrix_move_to_coo


    subroutine csr_sparse_matrix_move_from_coo(this, from)
    !-----------------------------------------------------------------
    !< Move from (COO) -> this (CSR)
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout) :: this
        class(coo_sparse_matrix_t), intent(inout) :: from
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
        call move_alloc(from%ia,itmp)
        call move_alloc(from%ja,this%ja)
        call move_alloc(from%val,this%val)
        call memalloc(this%get_num_rows()+1, this%irp, __FILE__, __LINE__)
        if(nnz <= 0) then
            this%irp(:) = 1
            return      
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
        call this%set_state(from%get_state())
        call memfree(itmp, __FILE__, __LINE__)
        call from%free()
    end subroutine csr_sparse_matrix_move_from_coo


    subroutine csr_sparse_matrix_move_to_fmt(this, to)
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
    end subroutine csr_sparse_matrix_move_to_fmt


    subroutine csr_sparse_matrix_move_from_fmt(this, from)
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
    end subroutine csr_sparse_matrix_move_from_fmt



    subroutine csr_sparse_matrix_apply_body(op,x,y) 
    !-----------------------------------------------------------------
    !< Apply matrix vector product
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(in)    :: op
        class(vector_t),            intent(in)    :: x
        class(vector_t) ,           intent(inout) :: y 
    !-----------------------------------------------------------------
        call x%GuardTemp()
        select type(x)
            class is (serial_scalar_array_t)
                select type(y)
                    class is(serial_scalar_array_t)
                        if (op%get_symmetric_storage()) then
                            call matvec_symmetric_storage(            &
                                        num_rows = op%get_num_rows(), &
                                        num_cols = op%get_num_cols(), &
                                        irp      = op%irp,            &
                                        ja       = op%ja,             &
                                        val      = op%val,            &
                                        x        = x%b,               &
                                        y        = y%b )
                        else
                            call matvec(num_rows = op%get_num_rows(), &
                                        num_cols = op%get_num_cols(), &
                                        irp      = op%irp,            &
                                        ja       = op%ja,             &
                                        val      = op%val,            &
                                        x        = x%b,               &
                                        y        = y%b )
                    end if
                end select
        end select
        call x%CleanTemp()
    contains

        subroutine matvec(num_rows, num_cols, irp, ja, val, x, y)
        !-------------------------------------------------------------
        !< Sparse matrix vector product
        !-------------------------------------------------------------
            integer(ip), intent(in)  :: num_rows
            integer(ip), intent(in)  :: num_cols
            integer(ip), intent(in)  :: irp(num_rows+1)
            integer(ip), intent(in)  :: ja(irp(num_rows+1)-1)
            real(rp)   , intent(in)  :: val(irp(num_rows+1)-1)
            real(rp)   , intent(in)  :: x(num_cols)
            real(rp)   , intent(out) :: y(num_rows)
            integer(ip)              :: ir,ic, iz
        !-------------------------------------------------------------
            y = 0.0_rp
            do ir = 1, num_rows
               do iz = irp(ir), irp(ir+1)-1
                  ic   = ja(iz)
                  y(ir) = y(ir) + x(ic)*val(iz)
               end do ! iz
            end do ! ir
        end subroutine matvec


        subroutine matvec_symmetric_storage(num_rows, num_cols, irp, ja, val, x, y)
        !-------------------------------------------------------------
        !< Symmetric stored sparse matrix vector product
        !-------------------------------------------------------------
            integer(ip), intent(in)  :: num_rows
            integer(ip), intent(in)  :: num_cols
            integer(ip), intent(in)  :: irp(num_rows+1)
            integer(ip), intent(in)  :: ja(irp(num_rows+1)-1)
            real(rp)   , intent(in)  :: val(irp(num_rows+1)-1)
            real(rp)   , intent(in)  :: x(num_cols)
            real(rp)   , intent(out) :: y(num_rows)
            integer(ip)              :: ir,ic, iz
        !-------------------------------------------------------------
            assert(num_rows==num_cols)
            y = 0.0_rp
            do ir = 1, num_rows
                y(ir) = y(ir) + x(ja(irp(ir)))*val(irp(ir))
                do iz = irp(ir)+1, irp(ir+1)-1
                    ic = ja(iz)
                    y(ir) = y(ir) + x(ic)*val(iz)
                    y(ic) = y(ic) + x(ir)*val(iz)
                end do ! iz
            end do ! ir
        end subroutine matvec_symmetric_storage

    end subroutine csr_sparse_matrix_apply_body


    subroutine csr_sparse_matrix_allocate_values_body(this, nz)
    !-----------------------------------------------------------------
    !< Allocate CSR values
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t), intent(inout)  :: this
        integer(ip), optional,      intent(in)     :: nz
    !-----------------------------------------------------------------
        check(.not. allocated(this%val))
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
                call this%insert(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
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
                call this%insert(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
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
                call this%insert(ia(i), ja(j), val(i+ioffset,j+joffset), 1, this%get_num_rows(), 1, this%get_num_cols())
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
                call this%insert(ia(i), ja(j), val(i+ioffset,j+joffset), 1, this%get_num_rows(), 1, this%get_num_cols())
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


    subroutine csr_sparse_matrix_split_2x2(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2
    !< A = [A_II A_IG]
    !<     [A_GI A_GG]
    !<
    !< this routine computes A_II, A_IG, A_GI and A_GG given the global 
    !< matrix A. Note that A_II, A_IG, A_GI and A_GG are all optional.
    !-----------------------------------------------------------------
        class(csr_sparse_matrix_t),                         intent(in)    :: this
        integer(ip),                                        intent(in)    :: num_row
        integer(ip),                                        intent(in)    :: num_col
        class(base_sparse_matrix_t), allocatable,           intent(inout) :: A_II
        class(base_sparse_matrix_t), allocatable,           intent(inout) :: A_IG
        class(base_sparse_matrix_t), allocatable, optional, intent(inout) :: A_GI
        class(base_sparse_matrix_t), allocatable,           intent(inout) :: A_GG
        integer(ip)                                                       :: i, j
        integer(ip)                                                       :: nz
        integer(ip)                                                       :: total_cols
        integer(ip)                                                       :: total_rows
    !-----------------------------------------------------------------
        assert(this%get_symmetric_storage() .eqv. A_II%get_symmetric_storage()) 
        assert(this%get_symmetric_storage() .eqv. A_GG%get_symmetric_storage())
        total_rows = this%get_num_rows()
        total_cols = this%get_num_cols()
        assert((A_II%get_num_rows()==num_row) .and. (A_II%get_num_cols()==num_col))
        assert((A_IG%get_num_rows()==num_row) .and. (A_IG%get_num_cols()==total_cols-num_col))
        assert((A_GG%get_num_rows()==total_rows-num_row) .and. (A_GG%get_num_cols()==total_cols-num_col))
        if(present(A_GI)) then
            assert((A_GI%get_num_rows()==this%get_num_rows()-num_row) .and. (A_GI%get_num_cols()==num_col))
        endif

        if(this%is_symbolic()) then
            do i=1, num_row
                do j=this%irp(i), this%irp(i+1)-1
                    if(this%ja(j)<=num_col) then
                        call A_II%insert(ia  = i,            &
                                         ja  = this%ja(j),   &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    else
                        call A_IG%insert(ia  = i,                  &
                                         ja  = this%ja(j)-num_col, &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    endif
                enddo
            enddo
            do i=num_row+1, total_rows
                do j=this%irp(i), this%irp(i+1)-1
                    if(this%ja(j)<=num_col) then
                        if(present(A_GI)) call A_GI%insert(ia  = i-num_row,    &
                                                           ja  = this%ja(j),   &
                                                           imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    else
                        call A_GG%insert(ia  = i-num_row,          &
                                         ja  = this%ja(j)-num_col, &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    endif
                enddo
            enddo
        else
            do i=1, num_row
                do j=this%irp(i), this%irp(i+1)-1
                    if(this%ja(j)<=num_col) then
                        call A_II%insert(ia  = i,            &
                                         ja  = this%ja(j),   &
                                         val = this%val(j),  &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    else
                        call A_IG%insert(ia  = i,                  &
                                         ja  = this%ja(j)-num_col, &
                                         val = this%val(j),        &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    endif
                enddo
            enddo
            do i=num_row+1, total_rows
                do j=this%irp(i), this%irp(i+1)-1
                    if(this%ja(j)<=num_col) then
                        if(present(A_GI)) call A_GI%insert(ia  = i-num_row,    &
                                                           ja  = this%ja(j),   &
                                                           val = this%val(j),  &
                                                           imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    else
                        call A_GG%insert(ia  = i-num_row,          &
                                         ja  = this%ja(j)-num_col, &
                                         val = this%val(j),        &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    endif
                enddo
            enddo
        endif
    end subroutine csr_sparse_matrix_split_2x2


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


end module csr_sparse_matrix_names

