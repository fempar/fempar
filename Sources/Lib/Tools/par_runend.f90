subroutine runend
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif 
    integer :: code, info, ierror
    logical :: initialized_mpi, finalized_mpi
    code = -1
    call mpi_initialized(initialized_mpi, ierror)
    call mpi_finalized(finalized_mpi, ierror)
#ifdef __GFORTRAN__
    call backtrace()
#endif 
    if(initialized_mpi .and. .not. finalized_mpi) then
        call mpi_abort(mpi_comm_world,code,info)
    else
        stop -1
    endif
end subroutine runend
