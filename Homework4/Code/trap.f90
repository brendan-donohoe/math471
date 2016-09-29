subroutine trap(upper, lower, u, jac, n, val)
  implicit none
  integer, intent(in) :: n
  real(kind = 8), intent(in) :: upper, lower, u(0:n), jac(0:n)
  real(kind = 8), intent(out) :: val
  real(kind = 8) :: x, h, fsum, flower, fupper
  !real(kind = 8), parameter :: pi = 4 * atan(1.d0)
  !real(kind = 8), external :: fun
  integer :: j

  fsum = 0.d0
  h = (upper - lower) / dble(n)
  do j = 1, n - 1
    fsum = fsum + u(lower + j * h)*jac(lower + j * h)
  end do
  flower = u(lower)
  fupper = u(upper)
  val = h * (((flower + fupper) / 2.d0) + fsum)
end subroutine trap