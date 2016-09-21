program gaussQuad
	implicit none
	double precision :: Integral_value,pi
	double precision, dimension(:), allocatable :: x, f, w
	integer :: n,s
	pi = acos(dble(-1.d0))
    do s = 1,2
        do n = 2,1000

            allocate(x(0:n),f(0:n),w(0:n))
            call lglnodes(x,w,n)
            f = exp(cos(x*(pi**dble(s))))
            Integral_value = sum(f*w)
            deallocate(x,f,w)
            print *, s, n, Integral_value
        end do
    end do
	
end program gaussQuad
