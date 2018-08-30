program par_test_h_adaptive_poisson
  use fempar_names
  use par_test_h_adaptive_poisson_driver_names
  implicit none
  integer(ip) :: i
  type(par_test_h_adaptive_poisson_fe_driver_t), save :: test_driver 
  type(mpi_context_t) :: world_context

  call fempar_init()  
  call test_driver%parse_command_line_parameters()
  ! It is possible to use an abstract context and allocate it
  ! at execution time after reading parameters.
  call world_context%create()
  call test_driver%setup_environment(world_context)
  call test_driver%setup_timers()
  do i = 1,1
    call test_driver%run_simulation()
    call test_driver%report_timers()
  end do
  call test_driver%free_timers()
  call test_driver%free_command_line_parameters()
  call test_driver%free_environment()
  call world_context%free(finalize=.true.)
  call fempar_finalize()
  
end program par_test_h_adaptive_poisson
