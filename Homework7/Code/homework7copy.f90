! In order to run the program with the specified number of processes, enter:
!
! mpirun -np (number of processes) homework7.x
!
! Once the program is run, each process gets its own copy of the program and runs through it.  Each is identifiable via its process id, numbered from 1 to (number of processes) - 1

program homework7
  use mpi
  use xycoord
  implicit none
  integer :: ierr, nprocs, myid ! ierr: error code for a hypothetical error we're never going to check.  nprocs: number of processes.  myid: the id of the process running the code.
  integer :: status(MPI_STATUS_SIZE) ! status: can pretty much ignore.
  
  integer, parameter :: nx = 91 ! nx: the total number of gridpoints in the x-direction.
  integer, parameter :: ny = 101 ! ny: the total number of gridpoints in the y-direction.

  integer :: nr,ns,nrl,nsl,i,j

  integer :: px,ix_off ! what is my x-proc index and what is my offset
  integer :: p_left,p_right,px_max
  integer :: p_down,p_up,py_max
  integer :: nxl,nyl !local size
  integer :: remx,remy
  real(kind = 8) :: hr,hs
  real(kind = 8), dimension(:),  allocatable :: r,s
  real(kind = 8), dimension(:,:),allocatable :: x,y,u,up,um,ur,us
  real(kind = 8), dimension(:,:),allocatable :: xr,xs,yr,ys,rx,ry,sx,sy,jac,a
  integer :: int_sum
  CHARACTER(7) :: charit

  call mpi_init(ierr) ! Initialize MPI stuff.
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr) ! Store the total number of processes in nprocs.
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr) ! Store the id of this process in myid.

  nr = nx ! Total number of gridpoints in the x-direction for the rectangular domain.
  ns = ny ! Total number of gridpoints in the y-direction for the rectangular domain.
  hr = 2.d0/dble(nr-1) ! Set our stepsize in the x-direction along the rectangular domain.
  hs = 2.d0/dble(ns-1) ! Same with y.

  ! Label the processes from 1 to px_max
  ! You will have to improve on this and split the
  ! grid in px_max*py_max = nprocs
  ! in such a way that communication is minimized
  px = myid + 1 ! The "index" of the current processor.
  px_max = nprocs ! The maximum index of all the processors.

  ! px = ???
  ! py = ???

  ! Split up the grid in the "x-direction"
  nxl = nx/px_max ! Assign the number of gridpoints this process has in the x-direction.
  remx = nx-nxl*px_max ! remx = nx % px_max
  if (px .le. remx) then ! If we fail to evenly divide the domain among all processes without having some leftover and this is one of the processes we give an extra gridpoint...
     nxl = nxl + 1 ! Give one extra gridpoint to this process.
     ix_off = (px-1)*nxl ! Assign the offset we'll add to in loop iterations from 1 to nxl in order to get this process's location on the main grid.  This, if this process does have an extra point.
  else
     ix_off = (remx)*(nxl+1) + (px-(remx+1))*nxl ! This, if it doesn't.
  end if
  ! write(*,*) 'Crazy arithmetic! : ',px,nxl,ix_off
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! Wait for all other processes to finish.
  call MPI_Reduce(nxl,int_sum,1,&
       MPI_INTEGER,MPI_SUM,&
       0,MPI_COMM_WORLD,ierr) ! Add up everyone's number of gridpoints in the x-direction and store it in int_sum (think this is just for debugging).
  if(myid == 0) then
     if (nx .ne. int_sum) then ! If the sum of everyone's number of x-gridpoints is not equal to the total number of gridpoints, stuff done went wrong!
        write(*,*) 'Something is wrong with the number of points in x-direction: ',&
             nx,int_sum
     end if
  end if
  nrl = nxl ! Number of r gridpoints we have on the square domain.
  ! You will also split the grid in the y-direction
  nyl = ny ! Number of y-gridpoints we're responsible for on the square domain.
  nsl = ns ! Number of s-gridpoints we're responsible for on the square domain.

  ! who are my neighbours? 
  ! NOTE THAT WE ARE MAPPING TO THE 1D INDEXING HERE!!!!
  p_left  = px-1 - 1 ! ID of the process in control of the domain to our left.
  p_right = px+1 - 1 ! ID of the process in control of the domain to our right.
  if (px .eq. px_max) p_right = MPI_PROC_NULL ! We have control of the rightmost part of the grid, we don't have a neighbor to our right.
  if (px .eq. 1) p_left = MPI_PROC_NULL ! We have control of the leftmost part of the grid, we don't have a neighbor to our left.

  ! Allocate memory for the various arrays
  allocate(r(1:nrl),s(1:nsl),u(0:nrl+1,0:nsl+1),up(0:nrl+1,0:nsl+1),um(0:nrl+1,0:nsl+1),&
       ur(0:nrl+1,0:nsl+1),us(0:nrl+1,0:nsl+1))
  allocate(x(0:nrl+1,0:nsl+1),y(0:nrl+1,0:nsl+1))
  allocate(xr(0:nrl+1,0:nsl+1),yr(0:nrl+1,0:nsl+1),xs(0:nrl+1,0:nsl+1),ys(0:nrl+1,0:nsl+1)&
       ,rx(0:nrl+1,0:nsl+1),ry(0:nrl+1,0:nsl+1),sx(0:nrl+1,0:nsl+1),sy(0:nrl+1,0:nsl+1)&
       ,jac(0:nrl+1,0:nsl+1),a(0:nrl+1,0:nsl+1))
  
  ! Give coordinates to each process.

  ! Looks like this is really similar to homework 4.
  do i = 1,nrl
     r(i) = -1.d0 + dble(i-1+ix_off)*hr ! Create a mapping of gridpoint indices to x-coordinates along the square domain for this particular process.
  end do
  do i = 1,nsl
     s(i) = -1.d0 + dble(i-1)*hs ! Same, but for y-coordinates.
  end do

  do j = 1,nsl
     do i = 1,nrl
        x(i,j) = x_coord(r(i),s(j)) ! Create a mapping of gridpoint indices to x-coordinates.
        y(i,j) = y_coord(r(i),s(j)) ! Create a mapping of gridpoint indices to y-coordinates.
     end do
  end do

  ! At a boundary we extrapolate
  x(1:nrl,0) = 2.d0*x(1:nrl,1)-x(1:nrl,2)
  y(1:nrl,0) = 2.d0*y(1:nrl,1)-y(1:nrl,2)
  x(1:nrl,nsl+1) = 2.d0*x(1:nrl,nsl)-x(1:nrl,nsl-1)
  y(1:nrl,nsl+1) = 2.d0*y(1:nrl,nsl)-y(1:nrl,nsl-1)

  ! Only proc 1 and px_max have boundaries
  if(px .eq. 1) then
     x(0,:) = 2.d0*x(1,:)-x(2,:)
     y(0,:) = 2.d0*y(1,:)-y(2,:)
  end if
  if(px .eq. px_max) then
     x(nrl+1,:) = 2.d0*x(nrl,:)-x(nrl-1,:)
     y(nrl+1,:) = 2.d0*y(nrl,:)-y(nrl-1,:)
  end if
  ! send to the left recieve from the right
  call MPI_Sendrecv(x(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,123,&
       x(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,123,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(y(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,124,&
       y(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,124,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(x(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,125,&
       x(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,125,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(y(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,126,&
       y(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,126,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! Communicate between the processes and send data between process domains.

  ! Now compute the metric
  xr(1:nrl,:) = (x(2:nrl+1,:)-x(0:nrl-1,:))/(2.d0*hr)
  yr(1:nrl,:) = (y(2:nrl+1,:)-y(0:nrl-1,:))/(2.d0*hr)
  xs(:,1:nsl) = (x(:,2:nsl+1)-x(:,0:nsl-1))/(2.d0*hs)
  ys(:,1:nsl) = (y(:,2:nsl+1)-y(:,0:nsl-1))/(2.d0*hs)

  ! At a boundary we extrapolate
  xr(1:nrl,0) = 2.d0*xr(1:nrl,1)-xr(1:nrl,2)
  xs(1:nrl,0) = 2.d0*xs(1:nrl,1)-xs(1:nrl,2)
  yr(1:nrl,0) = 2.d0*yr(1:nrl,1)-yr(1:nrl,2)
  ys(1:nrl,0) = 2.d0*ys(1:nrl,1)-ys(1:nrl,2)
  xr(1:nrl,nsl+1) = 2.d0*xr(1:nrl,nsl)-xr(1:nrl,nsl-1)
  xs(1:nrl,nsl+1) = 2.d0*xs(1:nrl,nsl)-xs(1:nrl,nsl-1)
  yr(1:nrl,nsl+1) = 2.d0*yr(1:nrl,nsl)-yr(1:nrl,nsl-1)
  ys(1:nrl,nsl+1) = 2.d0*ys(1:nrl,nsl)-ys(1:nrl,nsl-1)

  ! Only proc 1 and px_max have boundaries
  if(px .eq. 1) then
     xr(0,:) = 2.d0*xr(1,:)-xr(2,:)
     xs(0,:) = 2.d0*xs(1,:)-xs(2,:)
     yr(0,:) = 2.d0*yr(1,:)-yr(2,:)
     ys(0,:) = 2.d0*ys(1,:)-ys(2,:)
  end if
  if(px .eq. px_max) then
     xr(nrl+1,:) = 2.d0*xr(nrl,:)-xr(nrl-1,:)
     xs(nrl+1,:) = 2.d0*xs(nrl,:)-xs(nrl-1,:)
     yr(nrl+1,:) = 2.d0*yr(nrl,:)-yr(nrl-1,:)
     ys(nrl+1,:) = 2.d0*ys(nrl,:)-ys(nrl-1,:)
  end if

  ! send to the left recieve from the right
  call MPI_Sendrecv(xr(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,223,&
       xr(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,223,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(yr(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,224,&
       yr(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,224,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(xr(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,225,&
       xr(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,225,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(yr(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,226,&
       yr(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,226,MPI_COMM_WORLD,status,ierr)
  ! send to the left recieve from the right
  call MPI_Sendrecv(xs(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,323,&
       xs(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,323,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ys(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,324,&
       ys(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,324,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(xs(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,325,&
       xs(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,325,MPI_COMM_WORLD,status,ierr)
  call MPI_Sendrecv(ys(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,326,&
       ys(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,326,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  jac = xr*ys-xs*yr
  
  rx =  ys/jac
  ry = -xs/jac
  sx = -yr/jac
  sy =  xr/jac

  ! send to the left recieve from the right
  call MPI_Sendrecv(u(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,123,&
       u(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,123,MPI_COMM_WORLD,status,ierr)
  ! send to the right recieve from the left
  call MPI_Sendrecv(u(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,125,&
       u(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,125,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  u =  x + y 
  ! Now compute the derivative in the interior
  ur(1:nrl,1:nsl) = (u(2:nrl+1,1:nsl)-u(0:nrl-1,1:nsl))/(2.d0*hr)
  us(1:nrl,1:nsl) = (u(1:nrl,2:nsl+1)-u(1:nrl,0:nsl-1))/(2.d0*hs)
  up = rx*ur+sx*us
  WRITE(charit,"(I7.7)") myid
  call  printdble2d(x(1:nrl,1:nsl),nrl,nsl,'x'//charit//'.txt')
  call  printdble2d(y(1:nrl,1:nsl),nrl,nsl,'y'//charit//'.txt')
  call  printdble2d(up(1:nrl,1:nsl),nrl,nsl,'up'//charit//'.txt')
  
  
  ! Test Boundary Conditions
  ! Test Forcing
  ! Test Computation of right hand side 

  ! Finally do some..
  ! Time-stepping
  ! Output
  ! 

  ! Deallocate
  call mpi_finalize(ierr)
end program homework7

subroutine printdble2d(u,nx,ny,str)
  implicit none
  integer, intent(in) :: nx,ny
  real(kind = 8), intent(in) :: u(nx,ny)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do j=1,ny,1
     do i=1,nx,1
        write(2,fmt='(E24.16)',advance='no') u(i,j)
     end do
     write(2,'()')
  end do
  close(2)
end subroutine printdble2d
