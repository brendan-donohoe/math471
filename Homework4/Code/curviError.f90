subroutine curviError(ux, uy, uxEx, uyEx, jac, nr, ns, error)
	implicit none
	integer, intent(in) :: nr,ns
	real(kind = 8), intent(in) :: ux(0:nr,0:ns), uy(0:nr,0:ns)
	real(kind = 8), intent(in) :: uxEx(0:nr,0:ns), uxEx(0:nr,0:ns)
	real(kind = 8), intent(in) :: jac(0:nr,0:ns)
	real(kind = 8), intent(out) :: error
	real(kind = 8), dimension(:,:), allocatable :: funProduct
	real(kind = 8) :: val
	integer :: i,j

	allocate(funProduct(0:nr,0:ns))

	do j = 0,ns
		do i = 0,nr
			funProduct(i,j) = (ux(i,j) + uy(i,j) - uxEx(i,j) - uyEx(i,j))**2
			funProduct(i,j) = funProduct(i,j)*jac(i,j)
		end do
	end do

	error = 0.d0
	do i = 0,nr
		call trap(1.d0, -1.d0, funProduct(i,0:ns), ns, val)
		error = error + val
	end do
	error = (error * hr)**(0.5d0)
end subroutine curviError
