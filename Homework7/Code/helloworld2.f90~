! Basic program to demonstrate basic MPI functionality.  Each process sends an
! array containing the process id multiplied by 10 and process id multiplied by
! 20 to the next process.  Once finished, each process waits until everyone is
! done before leaving.

program helloworld
  use mpi
  implicit none
  integer :: ierr, nprocs, myid
  integer :: sendpid, recvpid
  integer, dimension(:), allocatable :: senddata, recvdata
  integer :: status(MPI_STATUS_SIZE)
  integer :: absmod

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

  allocate(senddata(1:2), recvdata(1:2))
  sendpid = absmod((myid + 1), nprocs)
  recvpid = absmod((myid - 1), nprocs)
  senddata(1) = myid * 10
  senddata(2) = myid * 20

  print *, "P", myid, ": Sent ", senddata, " to ", sendpid
  ! mpi_sendrecv(senddata, sendcount, sendtype, destid, sendtag, recvdata, recvcount, recvtype, recvid, recvtag, comm, status, ierr)
  ! Pretty much the send and recv methods in one method spliced at the end of tag from mpi_send and combined into one method.
  call mpi_sendrecv(senddata, 2, MPI_INTEGER, sendpid, 1, recvdata, 2, MPI_INTEGER, recvpid, 1, MPI_COMM_WORLD, status, ierr)
  print *, "P", myid, ": Received ", recvdata, " from ", recvpid
  
  ! mpi_barrier(comm, ierr)
  ! Blocks the process until all other processes reach this point.
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  print *, "Process ", myid, " finishes!"

  call mpi_finalize(ierr)
end program helloworld

function absmod(a, b) result(n)
  integer, intent(in) :: a, b
  integer :: n
  n = mod(a, b)
  if (n < 0) then
    n = n + b
  endif
end function absmod
