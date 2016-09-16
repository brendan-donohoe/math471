program trap
  implicit none
  real(kind = 8) :: x, h, fsum, flower, fupper, val
  real(kind = 8), parameter :: upper = 1.d0, lower = -1.d0, pi = 4 * atan(1.d0)
  real(kind = 8), external :: fun
  integer :: n, i, j
  integer, parameter :: nmax = 10000

  do i = 1, 2
    do n = 2, nmax
      fsum = 0.d0
      h = (upper - lower) / n
      do j = 1, n - 1
        fsum = fsum + fun(lower + j * h, pi**i)
      end do
      flower = fun(lower, pi**i)
      fupper = fun(upper, pi**i)
      val = h * (((flower + fupper) / 2.d0) + fsum)
      print *, i, n, val
    end do
  end do
end program trap

real(kind = 8) function fun(x, k)
  real(kind = 8), intent(in) :: x, k
  fun1 = exp(cos(k*x))
end function fun
