subroutine runend
#ifdef MPI_MOD
  use mpi
#endif
  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif	
  integer :: code, info 
  code = -1
  !write(output_unit,'(a)') 'Parallel !'
  call mpi_abort(mpi_comm_world,code,info)
end subroutine runend
