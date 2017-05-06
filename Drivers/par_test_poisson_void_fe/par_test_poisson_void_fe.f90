program par_test_poisson_void_fe
  use fempar_names
  use par_test_poisson_void_fe_driver_names
  !$ use omp_lib
  implicit none
  type(par_test_poisson_void_fe_fe_driver_t), save :: test_driver 
  !$OMP THREADPRIVATE(test_driver)

  !call sleep(20)
  !$OMP PARALLEL 
  !$ write(*,*) 'Begining with',omp_get_num_threads(),'threads'
  !$OMP BARRIER
  call fempar_init()  
  call test_driver%run_simulation()
  call fempar_finalize()
  !$OMP END PARALLEL   
end program par_test_poisson_void_fe
