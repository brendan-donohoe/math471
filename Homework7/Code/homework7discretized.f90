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
  
  integer, parameter :: nx = 91 ! nx: the total number of gridpoints in the x-direction. 91
  integer, parameter :: ny = 101 ! ny: the total number of gridpoints in the y-direction. 101
  real(kind = 8), parameter :: dt = 0.0001 ! dt: The stepsize in time used in time discretization.
  integer, parameter :: nsteps = 100 !Check at t = 0.8: Prev value: ~-15 compared to -0.72
  integer, parameter :: manutype = 0

  integer :: nr,ns,nrl,nsl,i,j,it

  integer :: px,ix_off ! what is my x-proc index and what is my offset
  integer :: p_left,p_right,px_max
  integer :: p_down,p_up,py_max
  integer :: nxl,nyl !local size
  integer :: remx,remy
  real(kind = 8) :: hr,hs
  real(kind = 8), dimension(:),  allocatable :: r,s
  real(kind = 8), dimension(:,:),allocatable :: x,y,u,up,um,ur,us
  real(kind = 8), dimension(:,:),allocatable :: xr,xs,yr,ys,rx,ry,sx,sy,jac,a
  real(kind = 8), dimension(:,:),allocatable :: Force,Lap,rhside,utt
  real(kind = 8) :: t
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

  ! Placeholders.
  p_up = MPI_PROC_NULL
  p_down = MPI_PROC_NULL

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

  !u =  x + y 
  ! Now compute the derivative in the interior
  !ur(1:nrl,1:nsl) = (u(2:nrl+1,1:nsl)-u(0:nrl-1,1:nsl))/(2.d0*hr)
  !us(1:nrl,1:nsl) = (u(1:nrl,2:nsl+1)-u(1:nrl,0:nsl-1))/(2.d0*hs)
  !up = rx*ur+sx*us
  !WRITE(charit,"(I7.7)") myid
  !call  printdble2d(x(1:nrl,1:nsl),nrl,nsl,'x'//charit//'.txt')
  !call  printdble2d(y(1:nrl,1:nsl),nrl,nsl,'y'//charit//'.txt')
  !call  printdble2d(up(1:nrl,1:nsl),nrl,nsl,'up'//charit//'.txt')
  
  ! Set u and um using the given initial conditions.
  call set_initial(um,u,nxl,nyl,x,y,dt,manutype)
  write(*,*) "Approx: u(", x(10,10), ",", y(10,10), ",0)=", um(10,10)
  write(*,*) "Actual: u(", x(10,10), ",", y(10,10), ",0)=", sin(2*0+acos(-1.d0)-x(10,10))*sin(y(10,10))
  !write(*,*) "Approx: u(", x(10,10), ",", y(10,10), ",", dt, ")=", u(10,10)
  !write(*,*) "Actual: u(", x(10,10), ",", y(10,10), ",", dt, ")=", sin(2*dt+acos(-1.d0)-x(10,10))sin(y(10,10))
  ! Set up a.
  call set_a(a,nxl,nyl,x,y,manutype)
  allocate(Force(0:nxl+1,0:nyl+1),rhside(0:nxl+1,0:nyl+1),utt(0:nxl+1,0:nyl+1),Lap(0:nxl+1,0:nyl+1))

  do it = 1,nsteps
    t = dt*dble(it)
    write(*,*) "Approx: u(", x(10,10), ",", y(10,10), ",", t, ")=", u(10,10)
    write(*,*) "Actual: u(", x(10,10), ",", y(10,10), ",", t, ")=", sin(2*t+acos(-1.d0)-x(10,10))*sin(y(10,10))
    ! Set the boundary conditions for u at the current time t.
    call set_bc(u,nxl,nyl,x,y,t,manutype,p_left,p_right,p_up,p_down)

    ! Send our left u points to the process to our left and receive right u
    ! points from the process to our right.
    call MPI_Sendrecv(u(1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,0,&
         u(nrl+1,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,0,&
         MPI_COMM_WORLD,status,ierr)
    ! Send right points to our right, receive left points from our left.
    call MPI_Sendrecv(u(nrl,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_right,0,&
         u(0,0:nsl+1),nsl+2,MPI_DOUBLE_PRECISION,p_left,0,&
         MPI_COMM_WORLD,status,ierr)
    ! Send lower points to below, receive upper points from above.
    call MPI_Sendrecv(u(0:nrl+1,1),nrl+2,MPI_DOUBLE_PRECISION,p_down,0,&
         u(0:nrl+1,nsl+1),nrl+2,MPI_DOUBLE_PRECISION,p_up,0,&
         MPI_COMM_WORLD,status,ierr)
    ! Send upper points to above, receive lower points from below.
    call MPI_Sendrecv(u(0:nrl+1,nsl),nrl+2,MPI_DOUBLE_PRECISION,p_up,0,&
         u(0:nrl+1,0),nrl+2,MPI_DOUBLE_PRECISION,p_down,0,&
         MPI_COMM_WORLD,status,ierr)

    ! Compute the forcing function for the current states of x and y and t.
    call compute_forcing(Force,nxl,nyl,x,y,t,manutype)
    ! Approximate the Laplacian using space discretization.
    call compute_lap(Lap,u,a,jac,rx,sx,ry,sy,nxl,nyl,hr,hs,t,x,y) !TODO: Delete last t parameter.
    !write(*,*) "Approximation: ", Lap(10,10)
    !rhside = Lap + Force
    call get_utt(rhside,nxl,nyl,x,y,t,manutype)

    !write(*,*) "rhside err: ", utt(10,10) - (Lap(10,10) + Force(10,10))
    ! Approximate the value of u at the next time step.
    up(1:nxl,1:nyl) = 2.d0*u(1:nxl,1:nyl)-um(1:nxl,1:nyl) + dt**2*rhside(1:nxl,1:nyl)
    !write(*,*) "u(", x(10,10), ",", y(10,10), ",", t+2*dt, ")=", up(10,10)

    ! Update the previous u to the current u, and the current u to the next u.
    um = u
    u = up
  end do

  call set_bc(u,nxl,nyl,x,y,t,manutype,p_left,p_right,p_up,p_down)
  write(*,*) "u(", x(10,10), ",", y(10,10), ",", t, ")=", u(10,10)

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

subroutine get_utt(utt,nxl,nyl,x,y,t,utt_type)
  implicit none
  integer, intent(in) :: nxl,nyl,utt_type
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: x,y
  real(kind = 8), intent(in) :: t
  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: utt
  
  real(kind = 8) :: omega,t0,kx,ky
  real(kind = 8), parameter :: pi = acos(-1.d0)

  if (utt_type == 0) then
    omega = 2
    t0 = pi / 2
    kx = 1
    ky = 1
    utt = -(omega**2)*sin(omega*(t+t0)-kx*x)*sin(ky*y)
  endif
end subroutine get_utt

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

subroutine set_a(a,nxl,nyl,x,y,a_type)
  implicit none
  integer, intent(in) :: nxl,nyl,a_type
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: x,y
  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: a
  
  if (a_type == 0) then
    a = x*y 
  endif
  !write(*,*) "x(10,10)=", x(10,10)
  !write(*,*) "y(10,10)=", y(10,10)
  !write(*,*) "a(10,10)=", a(10,10)
end subroutine set_a

subroutine set_initial(um,u,nxl,nyl,x,y,dt,ic_type)
  implicit none
  integer, intent(in) :: nxl,nyl,ic_type
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: x,y
  real(kind = 8), intent(in) :: dt
  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: um,u

  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: umt
  real(kind = 8) :: omega,t0,kx,ky
  real(kind = 8), parameter :: pi = acos(-1.d0)

  if (ic_type == 0) then ! Manufactured trigonometric solution.
    omega = 2
    t0 = pi / 2
    kx = 1
    ky = 1
    um = sin(omega*t0-kx*x)*sin(ky*y)
    umt = omega*cos(omega*t0-kx*x)*sin(ky*y)
    !u = sin(omega*(dt+t0)-kx*x)*sin(ky*y) !Exact value of u(x,y,dt) used to test convergence.
  endif
  u = um + umt*dt ! Euler method approximation of u at t=dt.
end subroutine set_initial

subroutine compute_forcing(Force,nxl,nyl,x,y,t,force_type)
  implicit none
  integer, intent(in) :: nxl,nyl,force_type
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: x,y
  real(kind = 8), intent(in) :: t
  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: Force

  real(kind = 8) :: omega,t0,kx,ky
  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: a,ax,ay,ux,uy,uxx,uyy,utt,lap
  real(kind = 8), parameter :: pi = acos(-1.d0)

  if (force_type == 0) then ! Manufactured trigonometric solution.
    ! u = sin(omega*(t+t_0)-k_x*x)*sin(k_y*y)
    ! a = x*y
    omega = 2
    t0 = pi / 2
    kx = 1
    ky = 1
    a = x*y
    ax = y
    ay = x
    ux = -kx*cos(omega*(t+t0)-kx*x)*sin(ky*y)
    uy = ky*sin(omega*(t+t0)-kx*x)*cos(ky*y)
    uxx = -(kx**2)*sin(omega*(t+t0)-kx*x)*sin(ky*y)
    uyy = -(ky**2)*sin(omega*(t+t0)-kx*x)*sin(ky*y)
    utt = -(omega**2)*sin(omega*(t+t0)-kx*x)*sin(ky*y)
  endif
  lap = ax*ux + a*uxx + ay*uy + a*uyy
  !write (*,*) "True Laplacian:", lap(10,10)
  Force = utt - lap
end subroutine compute_forcing

subroutine print_max_err(a,b,nx,ny)
  implicit none
  integer, intent(in) :: nx,ny
  real(kind = 8), dimension(0:nx+1,0:ny+1) :: a,b
  
  integer :: i,j,imax,jmax
  real(kind = 8) :: max_err,cur_err
  
  max_err = 0
  imax = 0
  jmax = 0
  do i = 0,nx+1
    do j = 0,ny+1
      cur_err = abs(a(i,j) - b(i,j))
      if (cur_err .gt. max_err) then
        max_err = cur_err
      endif
    end do
  end do
  write(*,*) "Maximum error: ", max_err
end subroutine

subroutine compute_lap(Lap,u,a,jac,rx,sx,ry,sy,nxl,nyl,hr,hs,t,x,y) !TODO: Delete last t parameter.
  implicit none
  integer, intent(in) :: nxl,nyl
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: u,a,jac,rx,sx,ry,sy
  real(kind = 8), intent(in) :: hr,hs
  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: Lap

  real(kind = 8), dimension(:,:,:), allocatable :: q,term
  integer :: i,j
  real(kind = 8), intent(in) :: t
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: x,y

  allocate(q(0:nxl+1,0:nyl+1,1:8),term(0:nxl+1,0:nyl+1,1:8))
  
  write(*,*) "At t=",t
  ! First, compute the different q functions.

  q(:,:,1) = jac*rx*a*rx
  q(:,:,2) = jac*rx*a*sx
  q(:,:,3) = jac*ry*a*ry
  q(:,:,4) = jac*ry*a*sy
  q(:,:,5) = jac*sx*a*rx
  q(:,:,6) = jac*sx*a*sx
  q(:,:,7) = jac*sy*a*ry
  q(:,:,8) = jac*sy*a*sy

  !write(*,*) "q1(10,10) approximation: ", q(10,10,1)
  !write(*,*) "q1(10,10) actual: ", x(10,10)*y(10,10)
  !write(*,*) "q8(10,9) approximation: ", q(10,9,8)
  !write(*,*) "q8(10,9) actual: ", x(10,9)*y(10,9)
  !write(*,*) "q8(10,10) approximation: ", q(10,10,8)
  !write(*,*) "q8(10,10) actual: ", x(10,10)*y(10,10)
  !write(*,*) "q8(10,11) approximation: ", q(10,11,8)
  !write(*,*) "q8(10,11) actual: ", x(10,11)*y(10,11)


  !write(*,*) "x(10,10): ", x(10,10)
  !write(*,*) "y(10,10): ", y(10,10)
  !write(*,*) "a(10,10): ", a(10,10)

  ! Next, use finite differences to approximate each term of the Laplacian.
  do i = 1,nxl
    do j = 1,nyl
      ! Term 1: (q1*ur) r
      term(i,j,1) = ((q(i+1,j,1) + q(i,j,1))*u(i+1,j) - (q(i+1,j,1) + 2*q(i,j,1) + q(i-1,j,1))*u(i,j) &
                    + (q(i,j,1) + q(i-1,j,1))*u(i-1,j)) / (2*hr**2)

      ! Term 2: (q2*us) r
      term(i,j,2) = ((u(i+1,j+1) - u(i+1,j-1))*q(i+1,j,2) - (u(i-1,j+1) - u(i-1,j-1))*q(i-1,j,2)) / (4*hs*hr)

      ! Term 3: (q3*ur) r
      term(i,j,3) = ((q(i+1,j,3) + q(i,j,3))*u(i+1,j) - (q(i+1,j,3) + 2*q(i,j,3) + q(i-1,j,3))*u(i,j) &
                    + (q(i,j,3) + q(i-1,j,3))*u(i-1,j)) / (2*hr**2)

      ! Term 4: (q4*us) r
      term(i,j,4) = ((u(i+1,j+1) - u(i+1,j-1))*q(i+1,j,4) - (u(i-1,j+1) - u(i-1,j-1))*q(i-1,j,4)) / (4*hs*hr)

      ! Term 5: (q5*ur) s
      term(i,j,5) = ((u(i+1,j+1) - u(i-1,j+1))*q(i,j+1,5) - (u(i+1,j-1) - u(i-1,j-1))*q(i,j-1,5)) / (4*hs*hr)

      ! Term 6: (q6*us) s
      term(i,j,6) = ((q(i,j+1,6) + q(i,j,6))*u(i,j+1) - (q(i,j+1,6) + 2*q(i,j,6) + q(i,j-1,6))*u(i,j) &
                    + (q(i,j,6) + q(i,j-1,6))*u(i,j-1)) / (2*hs**2)

      ! Term 7: (q7*ur) s
      term(i,j,7) = ((u(i+1,j+1) - u(i-1,j+1))*q(i,j+1,7) - (u(i+1,j-1) - u(i-1,j-1))*q(i,j-1,7)) / (4*hs*hr)

      ! Term 8: (q8*us) s
      term(i,j,8) = ((q(i,j+1,8) + q(i,j,8))*u(i,j+1) - (q(i,j+1,8) + 2*q(i,j,8) + q(i,j-1,8))*u(i,j) &
                    + (q(i,j,8) + q(i,j-1,8))*u(i,j-1)) / (2*hs**2)
    end do
  end do
  !write(*,*) "Finite difference approximations done."
  !write(*,*) "Approximation of term 1 at (i,j)=(10,10): ", term(10,10,1)
  !write(*,*) "Actual value of term 1 at (i,j)=(10,10): ", &
  !           y(10,10)*sin(y(10,10))*(cos(x(10,10)-2*t)-x(10,10)*sin(x(10,10)-2*t))
  !do i = 2,7
  !  write(*,*) "Approximation of term ", i, "at (i,j)=(10,10): ", term(10,10,i)
  !  write(*,*) "Actual value of of term ", i, " at (i,j)=(10,10): 0"
  !end do
  !write(*,*) "Approximation of term 8 at (i,j)=(10,10): ", term(10,10,8)
  !write(*,*) "Actual value of term 8 at (i,j)=(10,10): ", &
  !           x(10,10)*sin(x(10,10)-2*t)*(cos(y(10,10))-y(10,10)*sin(y(10,10)))
  Lap = term(:,:,1)
  do i = 2,8
    Lap = Lap + term(:,:,i)
  end do
  Lap = Lap / jac
  !write(*,*) "Approximation of Laplacian at (i,j)=(10,10): ", Lap(10,10)
  !write(*,*) "Actual value of Laplacian at (i,j)=(10,10): ", &
  !           y(10,10)*sin(y(10,10))*(cos(2*t-x(10,10))+x(10,10)*sin(2*t-x(10,10))) &
  !         + x(10,10)*sin(x(10,10)-2*t)*(cos(y(10,10))-y(10,10)*sin(y(10,10)))
  !write (*,*) "Value of sum of terms: ", (term(10,10,1) + term(10,10,2) + term(10,10,3) &
  !          + term(10,10,4) + term(10,10,5) + term(10,10,6) + term(10,10,7) + term(10,10,8))
  !  y*sin(y)*(cos(2*t-x)+x*sin(2*t-x)) + x*sin(x-2*t)*(cos(y)-y*sin(y))
  
  deallocate(q,term)
end subroutine compute_lap

subroutine set_bc(u,nxl,nyl,x,y,t,bc_type,p_left,p_right,p_up,p_down)
  use mpi
  implicit none
  integer, intent(in) :: nxl,nyl,bc_type
  integer, intent(in) :: p_left,p_right,p_up,p_down
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: x,y
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(inout) :: u
  real(kind = 8), intent(in) :: t

  real(kind = 8) :: omega,t0,kx,ky
  real(kind = 8), parameter :: pi = acos(-1.d0)
  
  ! First, let's make sure this process actually has some boundaries
  ! before doing any intensive work!
  if (p_left == MPI_PROC_NULL .or. p_right == MPI_PROC_NULL &
      .or. p_up == MPI_PROC_NULL .or. p_down == MPI_PROC_NULL) then
    if (bc_type == 0) then ! Manufactured trigonometric solution.
      omega = 2
      t0 = pi / 2
      kx = 1
      ky = 1
      if (p_left == MPI_PROC_NULL) then ! Leftmost boundary.
        u(1,:) = sin(omega*(t+t0)-kx*x(1,:))*sin(ky*y(1,:))
      endif
      if (p_right == MPI_PROC_NULL) then ! Rightmost boundary.
        u(nxl,:) = sin(omega*(t+t0)-kx*x(nxl,:))*sin(ky*y(nxl,:))
      endif
      if (p_down == MPI_PROC_NULL) then ! Lower boundary.
        u(:,1) = sin(omega*(t+t0)-kx*x(:,1))*sin(ky*y(:,1))
      endif
      if (p_up == MPI_PROC_NULL) then ! Upper boundary.
        u(:,nyl) = sin(omega*(t+t0)-kx*x(:,nyl))*sin(ky*y(:,nyl))
      endif
    endif
  endif
end subroutine set_bc
