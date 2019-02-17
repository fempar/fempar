program test_direct_solver

USE types_names
USE memor_names
USE fempar_names
USE sparse_matrix_names
USE direct_solver_names
USE direct_solver_parameters_names
USE direct_solver_creational_methods_dictionary_names
USE FPL
USE IR_Precision
use iso_c_binding


implicit none

# include "debug.i90"

    type(sparse_matrix_t)          :: sparse_matrix
    type(direct_solver_t)          :: direct_solver
    type(ParameterList_t)          :: parameter_list
    type(ParameterList_t), pointer :: direct_solver_params
    type(serial_scalar_array_t)    :: x
    type(serial_scalar_array_t)    :: y
    integer(ip), parameter         :: nrhs = 3
    integer(ip), parameter         :: n    = 4
    integer(ip), parameter         :: n_adapted = 3
    real(rp), allocatable          :: xx(:,:)
    real(rp), allocatable          :: yy(:,:)
    integer                        :: FPLError
    real(c_double)                 :: control_params(20) = 0
    integer                        :: iparm(64) = 0
    integer                        :: i, j, iters=5
    real(rp), parameter            :: exact_solution(n) = [ -0.3749999999999998E+00, &
                                                             0.3333333333333333E+00, &
                                                             -0.6250000000000001E+00, &
                                                             0.6874999999999999E+00 ]
    real(rp), parameter            :: tol = 1.0e-12                                                         
    real(rp), parameter            :: exact_solution_adapted_system(n_adapted) = [ -0.400000000000000, &
                                                                                    0.333333333333333, &
                                                                                    0.200000000000000 ]
    call meminit()
    ! ParameterList: initialize
    call FPL_Init()
    call the_direct_solver_creational_methods_dictionary%Init()
    call parameter_list%Init()

    call create_sparse_matrix()
    call create_arrays()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PARDISO MKL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    direct_solver_params => parameter_list%NewSubList(Key=pardiso_mkl)

    ! ParameterList: set parameters
    FPLError = 0
    FPLError = FPLError + direct_solver_params%set(key = dls_type_key,        value = pardiso_mkl)
    FPLError = FPLError + direct_solver_params%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_uss)
    FPLError = FPLError + direct_solver_params%set(key = pardiso_mkl_message_level, value = 0)
    FPLError = FPLError + direct_solver_params%set(key = pardiso_mkl_iparm,         value = iparm)
    check(FPLError == 0)

    do i=1, iters
#ifdef ENABLE_MKL
      ! Direct solver: create and set properties
      if(i==1) then
          call direct_solver%set_type_from_pl(direct_solver_params)
          call direct_solver%set_matrix(sparse_matrix)
      else
          call direct_solver%set_parameters_from_pl(direct_solver_params)
      endif
        
      call solve_and_check_solution_single_rhs(exact_solution)
      call solve_and_check_solution_multiple_rhs()

      if(i/=iters) then
          print*, ''
          print*, '!< =============================================='
          print*, '!< UPDATE MATRIX WITH SAME_NONZERO_PATTERN:', mod(i,2)==0
          print*, '!< =============================================='
      endif
        
      call scale_sparse_matrix_and_associated_arrays()

      call direct_solver%replace_matrix(sparse_matrix, same_nonzero_pattern=mod(i,2)==0)
      call direct_solver%free_in_stages(free_numerical_setup)
      call direct_solver%update_matrix(same_nonzero_pattern=mod(i,2)==0)

      if ( i == iters ) then
        ! Test update_after_remesh TBP of direct_solver_t
        ! 1. Destroy sparse matrix and associated arrays and 
        !    generate new ones of different size
        call create_adapted_sparse_matrix()
        call create_arrays()

        ! 2. Update direct solver to reflect that the sparse matrix has changed 
        !    (typically after AMR mesh adaptation, although not necessarily)
        call direct_solver%reallocate_after_remesh()
        call solve_and_check_solution_single_rhs(exact_solution_adapted_system)
      end if 
#endif
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! UMFPACK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    direct_solver_params => parameter_list%NewSubList(Key=UMFPACK)

    ! ParameterList: set parameters
    FPLError = 0
    FPLError = FPLError + direct_solver_params%set(key = dls_type_key,     value = UMFPACK)
    FPLError = FPLError + direct_solver_params%set(key = UMFPACK_CONTROL_PARAMS, value = control_params)
    check(FPLError == 0)

    do i=1, iters

#ifdef ENABLE_UMFPACK
        ! Direct solver: create and set properties
        if(i==1) then
            call direct_solver%set_type_from_pl(direct_solver_params)
            call direct_solver%set_matrix(sparse_matrix)
        else
            call direct_solver%set_parameters_from_pl(direct_solver_params)
        endif

        call solve_and_check_solution_single_rhs()
        call solve_and_check_solution_multiple_rhs()

        if(i/=iters) then
            print*, ''
            print*, '!< =============================================='
            print*, '!< UPDATE MATRIX WITH SAME_NONZERO_PATTERN:', mod(i,2)==0
            print*, '!< =============================================='
        endif
        
        call scale_sparse_matrix_and_associated_arrays()

        call direct_solver%replace_matrix(sparse_matrix, same_nonzero_pattern=mod(i,2)==0)
        call direct_solver%free_in_stages(free_numerical_setup)
        call direct_solver%update_matrix(same_nonzero_pattern=mod(i,2)==0)
#endif

    enddo

    ! Free
    call parameter_list%free()
    call direct_solver%free()
    call sparse_matrix%free()
    call x%free()
    call y%free()
    call memfree(xx, __FILE__, __LINE__)
    call memfree(yy, __FILE__, __LINE__)

    call memstatus()
    
contains 
   function nrm2_error ( exact, approx ) 
      implicit none
      real(rp), intent(in) :: exact(:)
      real(rp), intent(in) :: approx(:)
      real(rp) :: nrm2_error
      integer(ip) :: i 
      assert (size(exact) == size(approx))
      nrm2_error = 0.0 
      do i=1, size(exact)
        nrm2_error = nrm2_error + (exact(i)-approx(i))**2.0_rp
      end do 
      nrm2_error = sqrt(nrm2_error)
   end function nrm2_error

   subroutine create_sparse_matrix()
     implicit none
     call sparse_matrix%free()
     ! Sparse matrix: create
     call sparse_matrix%create(n, .false., .true., SPARSE_MATRIX_SIGN_UNKNOWN )
     call sparse_matrix%insert(nz=8,                    &
                               ia=[1,1,2,3,3,4,4,4],  &
                               ja=[1,4,2,3,4,1,3,4],  &
                               val=[1.,2.,3.,5.,6.,2.,6.,8.])
     call sparse_matrix%convert('CSR')
     call sparse_matrix%print(6)
     call sparse_matrix%set_sum_duplicates(.false.)
   end subroutine  create_sparse_matrix 

   subroutine create_arrays()
     implicit none
     ! Serial scalar array: create and initialize
     call x%free()
     call y%free()
     call x%create_and_allocate(sparse_matrix%get_num_rows())
     call x%init(1.0_rp)
     call y%create_and_allocate(sparse_matrix%get_num_cols())
     call x%print(6)
     if (allocated(xx)) then
       call memfree(xx, __FILE__, __LINE__)
     end if  
     call memalloc(sparse_matrix%get_num_rows(), nrhs, xx, __FILE__, __LINE__)
     xx = 1.0_rp
     if (allocated(yy)) then
       call memfree(yy, __FILE__, __LINE__)
     end if  
     call memalloc(sparse_matrix%get_num_cols(), nrhs, yy, __FILE__, __LINE__)
   end subroutine create_arrays

   subroutine create_adapted_sparse_matrix()
     implicit none
     call sparse_matrix%free()
     call sparse_matrix%create(n_adapted, .false., .true., SPARSE_MATRIX_SIGN_UNKNOWN )
     call sparse_matrix%insert(nz=4, &
                              ia=[1,1,2,3],  &
                              ja=[1,3,2,3],  &
                              val=[1.,7.,3.,5.])
     call sparse_matrix%convert('CSR')
     call sparse_matrix%print(6)
  end subroutine  create_adapted_sparse_matrix 
   
   subroutine scale_sparse_matrix_and_associated_arrays()
     implicit none
     assert(allocated(xx))
     ! Scale both sizes of the system by diagonal matrix D_ii = {i+1}
     ! The solution must remain the same
     call sparse_matrix%insert(nz=8,                    &
                              ia=[1,1,2,3,3,4,4,4],  &
                              ja=[1,4,2,3,4,1,3,4],  &
                              val=[1.,2.,3.,5.,6.,2.,6.,8.]*(i+1))
     call sparse_matrix%convert('CSR')
     call sparse_matrix%print(6)
     call x%init(real(i+1,rp))
     xx(:,:) = real(i+1,rp)
   end subroutine scale_sparse_matrix_and_associated_arrays

   subroutine solve_and_check_solution_multiple_rhs()
     implicit none
     call direct_solver%solve(xx,yy) 
     call direct_solver%log_info()
     do j=1, nrhs
        write(*,*) 'i=', j, '||x_i-x_i*||_2=', nrm2_error(exact_solution,yy(:,j)) 
        check ( nrm2_error(exact_solution,yy(:,j)) < tol )
     end do 
   end subroutine solve_and_check_solution_multiple_rhs 

   subroutine solve_and_check_solution_single_rhs(exact_solution)
     implicit none
     real(rp), intent(in) :: exact_solution(:)
     ! Direct solver: analisys, factorization and solve
     call direct_solver%solve(x,y)
     call direct_solver%log_info()
     write(*,*) '||x-x*||_2=', nrm2_error(exact_solution,y%get_entries()) 
     check ( nrm2_error(exact_solution,y%get_entries()) < tol )
   end subroutine solve_and_check_solution_single_rhs 

end program test_direct_solver
