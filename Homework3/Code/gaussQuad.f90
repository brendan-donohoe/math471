program gaussQuad
	implicit none
	double precision :: f,Integral_value,x,pi
	integer :: w,n,s
	w = 1
	pi = acos(dble(-1.d0))
    do s = 1,2
        do n = 2,1000

            allocate(x(0:n),f(0:n),w(0:n))
            call lglnodes(x,w,n)
            f = exp(cos(x*(pi**dble(s)))
            Integral_value = sum(f*w)
            deallocate(x,f,w)
        end do
    end do
	
end program gaussQuad