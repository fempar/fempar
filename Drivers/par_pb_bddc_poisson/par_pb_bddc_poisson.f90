program par_pb_bddc_poisson
  use fempar_names
  use par_pb_bddc_poisson_driver_names
  implicit none
  integer(ip) :: i
  type(par_pb_bddc_poisson_fe_driver_t)     :: test_driver  
  type(mpi_context_t) :: world_context

  call world_context%create()
  call fempar_init()  
  call test_driver%parse_command_line_parameters()
  call test_driver%setup_environment(world_context)
  do i = 1,1
     call test_driver%run_simulation()
  end do 
  call test_driver%free_command_line_parameters()
  call test_driver%free_environment()
  call fempar_finalize()
  call world_context%free(finalize=.true.)
end program par_pb_bddc_poisson
