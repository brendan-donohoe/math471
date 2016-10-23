program hwk4
  !$ use omp_lib
  use xycoord ! use the module xycoord to set the mapping
  implicit none

  integer :: nr,ns,i,j

  real(kind = 8) :: hr,hs,det,integral,val,error,tstart,tend
  real(kind = 8), dimension(:,:), allocatable :: jacproduct

  real(kind = 8), dimension(:), allocatable :: r,s
  real(kind = 8), dimension(:,:), allocatable :: u
  real(kind = 8), dimension(:,:), allocatable :: xc,yc
  real(kind = 8), dimension(:,:), allocatable :: xr,xs
  real(kind = 8), dimension(:,:), allocatable :: yr,ys
  real(kind = 8), dimension(:,:), allocatable :: jac

  !$ call OMP_set_num_threads(1)

  nr = 800
  ns = 800
  
  ! Allocate memory for the various arrays
  allocate(r(0:nr),s(0:ns),u(0:nr,0:ns))
  allocate(xc(0:nr,0:ns),yc(0:nr,0:ns))
  allocate(xr(0:nr,0:ns),xs(0:nr,0:ns),yr(0:nr,0:ns),ys(0:nr,0:ns))
  allocate(jac(0:nr,0:ns))
  allocate(jacproduct(0:nr,0:ns))
  
  tstart = omp_get_wtime()

  hr = 2.d0/dble(nr)
  hs = 2.d0/dble(ns)
  !$OMP PARALLEL DO PRIVATE(i)
  do i = 0,nr
     r(i) = -1.d0 + dble(i)*hr
  end do
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO PRIVATE(i)
  do i = 0,ns
     s(i) = -1.d0 + dble(i)*hs
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO PRIVATE(j,i)
  do j = 0,ns
     do i = 0,nr
        xc(i,j) = x_coord(r(i),s(j))
        yc(i,j) = y_coord(r(i),s(j))
     end do
  end do
  !$OMP END PARALLEL DO
  !call  printdble2d(xc,nr,ns,'x.txt')
  !call  printdble2d(yc,nr,ns,'y.txt')

  ! Differentiate x and y with respect to r
  !$OMP PARALLEL DO PRIVATE(i)
  do i = 0,ns
    call differentiate(xc(0:nr,i),xr(0:nr,i),hr,nr)
    call differentiate(yc(0:nr,i),yr(0:nr,i),hr,nr)
  end do
  !$OMP END PARALLEL DO

  ! Differentiate x and y with respect to s
  !$OMP PARALLEL DO PRIVATE(i)
  do i = 0,nr
    call differentiate(xc(i,0:ns),xs(i,0:ns),hs,ns)
    call differentiate(yc(i,0:ns),ys(i,0:ns),hs,ns)
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO PRIVATE(j,i)
  do j = 0,ns
     do i = 0,nr
      !u(i,j) = exp(xc(i,j)+yc(i,j))
      u(i,j) = 1.d0
     end do
  end do
  !$OMP END PARALLEL DO

  ! Calculate Jacobian (yay!)
  !$OMP PARALLEL DO PRIVATE(j,i)
  do j = 0, ns
    do i = 0, nr
      jac(i,j) = xr(i,j)*ys(i,j) - xs(i,j)*yr(i,j)
    end do
  end do
  !$OMP END PARALLEL DO

  ! Store the elements of the Jacobian multiplied by the elements of u.
  !$OMP PARALLEL DO PRIVATE(j,i)
  do j = 0,ns
    do i = 0,nr
      jacproduct(i,j) = u(i,j)*jac(i,j)
    end do
  end do
  !$OMP END PARALLEL DO

  ! Integrate (actually estimate using trapezoidal rule)!
  integral = 0.d0
  !$OMP PARALLEL DO PRIVATE(i,val) REDUCTION(+:integral)
  do i = 0, nr
    call trap(1.d0, -1.d0, jacproduct(i,0:ns), ns, val)
    integral = integral + val
  end do
  !$OMP END PARALLEL DO
  integral = integral*hr

  tend = omp_get_wtime()
  
  write(*,*) "integral value=", integral
  write(*,*) "time taken=", tend - tstart
end program hwk4
