program par_test_nsi
  use fempar_names
  use par_test_nsi_driver_names
  !$ use omp_lib
  implicit none
  integer(ip) :: i
  type(par_test_nsi_fe_driver_t), save :: test_driver 

#ifdef _OPENMP   
  type(mpi_omp_context_t) :: world_context
#else
  type(mpi_context_t) :: world_context
#endif
  
  !$OMP THREADPRIVATE(test_driver)

  !call sleep(20)
  !$OMP PARALLEL 
  !$ write(*,*) 'Begining with',omp_get_num_threads(),'threads'
  !$OMP BARRIER
  call world_context%create()
  call fempar_init()  
  call test_driver%parse_command_line_parameters()
  call test_driver%setup_environment(world_context)
  call test_driver%setup_timers()
  do i = 1,1
    call test_driver%run_simulation()
    !call test_driver%report_timers()
  end do
  !call test_driver%free_timers()
  call test_driver%free_command_line_parameters()
  call test_driver%free_environment()
  call fempar_finalize()
  call world_context%free(finalize=.true.)
  !$OMP END PARALLEL   
end program par_test_nsi
