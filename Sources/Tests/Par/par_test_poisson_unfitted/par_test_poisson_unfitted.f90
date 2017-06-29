program par_test_poisson_unfitted
  use fempar_names
  use par_test_poisson_unfitted_driver_names
  !$ use omp_lib
  implicit none
  integer(ip) :: i
  type(par_test_poisson_unfitted_fe_driver_t), save :: test_driver 
  !$OMP THREADPRIVATE(test_driver)

  !call sleep(20)
  !$OMP PARALLEL 
  !$ write(*,*) 'Begining with',omp_get_num_threads(),'threads'
  !$OMP BARRIER
  call fempar_init()  
  call test_driver%parse_command_line_parameters()
  call test_driver%setup_environment()
  call test_driver%setup_timers()
  do i = 1,test_driver%test_params%get_number_runs()
    call test_driver%run_simulation()
    call test_driver%report_timers()
  end do
  call test_driver%free_timers()
  call test_driver%free_command_line_parameters()
  call test_driver%free_environment()
  call fempar_finalize()
  !$OMP END PARALLEL   
end program par_test_poisson_unfitted
