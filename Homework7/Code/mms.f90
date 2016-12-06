module mms
  real(kind = 8), parameter :: pi = acos(-1.d0)
contains

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

subroutine set_a(a,nxl,nyl,x,y,a_type)
  implicit none
  integer, intent(in) :: nxl,nyl,a_type
  real(kind = 8), dimension(0:nxl+1,0:nyl+1), intent(in) :: x,y
  real(kind = 8), dimension(0:nxl+1,0:nyl+1) :: a
  
  if (a_type == 0) then
    a = x*y 
  endif
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
  Force = utt - lap
end subroutine compute_forcing

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
