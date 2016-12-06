! In order to run the program with the specified number of processes, enter:
!
! mpirun -np (number of processes) helloworld.x
!
! Basic program to demonstrate basic MPI functionality.  Each process sends its
! process id multiplied by 10 to the process id after it (so, 0 sends 0 to 1, 1
! sends 10 to 2, ..., 7 sends 70 to 0).

program helloworld
  use mpi
  implicit none
  integer :: ierr, nprocs, myid
  integer :: status(MPI_STATUS_SIZE) ! Can pretty much ignore this.
  integer :: absmod
  integer :: sendpid, recvpid, senddata, recvdata

  call mpi_init(ierr) ! Enable us to use MPI functionality.
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr) ! Get the total number of processes running and store it in nprocs.
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr) ! Get the particular process running this instance of the program's id.

  sendpid = absmod((myid + 1), nprocs)
  recvpid = absmod((myid - 1), nprocs)
  senddata = myid * 10

  print *, "P", myid, ": Sent ", senddata, " to ", sendpid

  ! mpi_send(data, count, type, destid, tag, comm, ierr)
  ! data: The data we're sending to the process corresponding to the "type" below.
  ! count: The number of elements we're sending, if data is an array.  In this case, the data is sent starting from index 1.
  ! type: the type of the data we're sending (handlenames as defined by MPI; such as: MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION).
  ! destid: The id of the process that will be receiving the data.
  ! tag: The "mailbox" of the process we'll be sending to which the process will check for new data.
  ! comm: The communicator, pretty much always MPI_COMM_WORLD.
  ! ierr: The error code returned if something went wrong with sending the data.
  call mpi_send(senddata, 1, MPI_INTEGER, sendpid, 1, MPI_COMM_WORLD, ierr)


  ! mpi_recv(data, count, type, sendid, tag, comm, status, ierr)
  ! data: The variable in which the received data is stored.
  ! count: The number of elements we're receiving, if data is an array.  In this case, the data is received starting from index 1.
  ! type: the type of the data we're receiving.
  ! destid: The id of the process that has sent the data.
  ! tag: The "mailbox" in which the process looks for information.
  ! comm: The communicator, still pretty much always MPI_COMM_WORLD.
  ! status: Variable in which technical stuff about the transfer is stored, can pretty much be ignored.
  ! ierr: The error code returned if something went wrong with sending the data.
  call mpi_recv(recvdata, 1, MPI_INTEGER, recvpid, 1, MPI_COMM_WORLD, status, ierr)
  print *, "P", myid, ": Received ", recvdata, " from ", recvpid

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
