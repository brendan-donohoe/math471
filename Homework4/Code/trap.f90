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
