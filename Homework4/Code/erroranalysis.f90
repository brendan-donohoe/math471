program erroranalysis
  implicit none
  integer :: k

  do k = 5, 400
    call computeerror(k, 2*k)
  end do
end program

subroutine computeerror(nr, ns)
  use xycoord ! use the module xycoord to set the mapping 
  implicit none

  integer, intent(in) :: nr, ns

  integer :: i,j,k,maxk

  real(kind = 8) :: hr,hs,det,val,error,jmax
  real(kind = 8), dimension(:,:), allocatable :: jacproduct

  real(kind = 8), dimension(:), allocatable :: r,s
  real(kind = 8), dimension(:,:), allocatable :: u,ur,us
  real(kind = 8), dimension(:,:), allocatable :: xc,yc
  real(kind = 8), dimension(:,:), allocatable :: xr,xs
  real(kind = 8), dimension(:,:), allocatable :: yr,ys
  real(kind = 8), dimension(:,:), allocatable :: rx,ry
  real(kind = 8), dimension(:,:), allocatable :: sx,sy
  real(kind = 8), dimension(:,:), allocatable :: ux,uy
  real(kind = 8), dimension(:,:), allocatable :: jac
  real(kind = 8), dimension(:,:), allocatable :: uxex, uyex

  ! Allocate memory for the various arrays
  allocate(r(0:nr),s(0:ns),u(0:nr,0:ns),ur(0:nr,0:ns),us(0:nr,0:ns))
  allocate(xc(0:nr,0:ns),yc(0:nr,0:ns))
  allocate(xr(0:nr,0:ns),xs(0:nr,0:ns),yr(0:nr,0:ns),ys(0:nr,0:ns))
  allocate(rx(0:nr,0:ns),ry(0:nr,0:ns),sx(0:nr,0:ns),sy(0:nr,0:ns))
  allocate(ux(0:nr,0:ns),uy(0:nr,0:ns))
  allocate(jac(0:nr,0:ns))
  allocate(uxex(0:nr,0:ns),uyex(0:nr,0:ns))
  allocate(jacproduct(0:nr,0:ns))

  hr = 2.d0/dble(nr)
  hs = 2.d0/dble(ns)
  do i = 0,nr
    r(i) = -1.d0 + dble(i)*hr
  end do
  do i = 0,ns
    s(i) = -1.d0 + dble(i)*hs
  end do
  
  do j = 0,ns
    do i = 0,nr
      xc(i,j) = x_coord(r(i),s(j))
      yc(i,j) = y_coord(r(i),s(j))
    end do
  end do

  ! Differentiate x and y with respect to r
  do i = 0,ns
    call differentiate(xc(0:nr,i),xr(0:nr,i),hr,nr)
    call differentiate(yc(0:nr,i),yr(0:nr,i),hr,nr)
  end do
  ! Differentiate x and y with respect to s
  do i = 0,nr
    call differentiate(xc(i,0:ns),xs(i,0:ns),hs,ns)
    call differentiate(yc(i,0:ns),ys(i,0:ns),hs,ns)
  end do

  ! Calculate the metric rx, ry, sx, sy
  do j = 0,ns
    do i = 0,nr
        det = xr(i,j)*ys(i,j) - yr(i,j)*xs(i,j)
        rx(i,j) = ys(i,j)/det
        ry(i,j) = -xs(i,j)/det
        sx(i,j) = -yr(i,j)/det
        sy(i,j) = xr(i,j)/det
    end do
  end do
  
  do j = 0,ns
    do i = 0,nr
      ! FUNCTION 0:
      u(i,j) = 1.d0
      uxex(i,j) = 0.d0
      uyex(i,j) = 0.d0
  
      ! FUNCTION 1:
      ! u(i,j) = sin(xc(i,j))*cos(yc(i,j))
      ! uxex(i,j) = cos(xc(i,j))*cos(yc(i,j))
      ! uyex(i,j) = -sin(xc(i,j))*sin(yc(i,j))

      ! FUNCTION 2:
      ! u(i,j) = exp(xc(i,j)+yc(i,j))
      ! uxex(i,j) = exp(xc(i,j)+yc(i,j))
      ! uyex(i,j) = exp(xc(i,j)+yc(i,j))

      ! FUNCTION 3:
      ! u(i,j) = xc(i,j)**2.d0 + yc(i,j)**2.d0
      ! uxex(i,j) = 2.d0*xc(i,j)
      ! uyex(i,j) = 2.d0*yc(i,j)
    end do
  end do
  
  ! Differentiate in the r-direction
  do i = 0,ns
    call differentiate(u(0:nr,i),ur(0:nr,i),hr,nr)
  end do
  ! Differentiate in the s-direction
  do i = 0,nr
    call differentiate(u(i,0:ns),us(i,0:ns),hs,ns)
  end do

  ! Calculate the approximate derivative of u with respect to x and y.
  do j = 0,ns
    do i = 0,nr
      ux(i,j) = rx(i,j)*ur(i,j) + sx(i,j)*us(i,j)
      uy(i,j) = ry(i,j)*ur(i,j) + sy(i,j)*us(i,j)
    end do
  end do

  jmax = xr(0,0)*ys(0,0) - xs(0,0)*yr(0,0)

  ! Calculate Jacobian and get its maximum value.
  do j = 0, ns
    do i = 0, nr
      jac(i,j) = xr(i,j)*ys(i,j) - xs(i,j)*yr(i,j)
      if (jac(i,j) > jmax) then
        jmax = jac(i,j)
      endif
    end do
  end do

  call curvierror(ux, uy, uxex, uyex, jac, nr, ns, error)
  print *, (hr*hs*jmax)**(1.d0/2), " ", error
end subroutine computeerror
