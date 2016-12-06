subroutine curvierror(u, v, jac, nrl, nsl, error)
  implicit none
  integer, intent(in) :: nrl,nsl
  real(kind = 8), intent(in) :: u(0:nrl+1,0:nsl+1), v(0:nrl+1,0:nsl+1)
  real(kind = 8), intent(in) :: jac(0:nrl+1,0:nsl+1)
  real(kind = 8), intent(out) :: error

  real(kind = 8), dimension(:,:), allocatable :: interior
  real(kind = 8) :: trapval
  integer :: i

  allocate(interior(1:nrl,1:nsl))

  interior = ((u(1:nrl,1:nsl) - v(1:nrl,1:nsl))**2)*jac(1:nrl,1:nsl)

  error = 0.d0
  do i = 1,nrl
    call trap(1.d0, -1.d0, interior(i,:), nsl, trapval)
    error = error + trapval
  end do
  error = (error * (2.d0/dble(nr)))**(0.5d0)
end subroutine curvierror

subroutine trap(upper, lower, f, n, val)
  implicit none
  integer, intent(in) :: n
  real(kind = 8), intent(in) :: upper, lower, f(0:n)
  real(kind = 8), intent(out) :: val

  integer :: i
  real(kind = 8) :: dx, fsum
  
  dx = (upper - lower) / n
  fsum = 0
  do i = 1,n-1
    fsum = fsum + f(i)
  end do

  val = (dx / 2.d0) * (f(0) + f(n) + 2.d0*fsum)
end subroutine trap
