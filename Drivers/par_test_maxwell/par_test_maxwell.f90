program par_test_maxwell
  use fempar_names
  use par_test_maxwell_driver_names
  implicit none
  integer(ip) :: i
  type(par_test_maxwell_fe_driver_t), save :: test_driver 

  call fempar_init()  
  call test_driver%parse_command_line_parameters()
  call test_driver%setup_environment()
  call test_driver%setup_timers()
  do i = 1,1
    call test_driver%run_simulation()
    call test_driver%report_timers()
  end do
  call test_driver%free_timers()
  call test_driver%free_command_line_parameters()
  call test_driver%free_environment()
  call fempar_finalize()
 
end program par_test_maxwell
