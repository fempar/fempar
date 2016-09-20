program par_test_poisson
  use fempar_names
  use par_test_poisson_driver_names
  implicit none
  type(par_test_poisson_fe_driver_t)     :: test_driver  
  call fempar_init()  
  call test_driver%run_simulation()
  call fempar_finalize()
end program par_test_poisson
