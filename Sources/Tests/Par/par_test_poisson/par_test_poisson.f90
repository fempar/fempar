program par_test_poisson
  use fempar_names
  use par_test_poisson_driver_names
  !$ use omp_lib
  implicit none
  type(par_test_poisson_fe_driver_t), save :: test_driver 
  !$OMP THREADPRIVATE(test_driver)

  call sleep(30)  
  !$OMP PARALLEL 
  !!!$OMP PARALLEL PRIVATE(test_driver) 
  !$ write(*,*) 'Begining with',omp_get_num_threads(),'threads'
  !$OMP BARRIER
  !call sleep(10)
  call fempar_init()  
  call test_driver%run_simulation()
  call fempar_finalize()
  !$OMP END PARALLEL   
end program par_test_poisson
