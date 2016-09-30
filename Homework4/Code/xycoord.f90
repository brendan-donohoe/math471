module xycoord
  real(kind = 8), parameter :: pi = acos(-1.d0)
  save
contains
  
  real(kind=8) function x_coord(r,s)
    implicit none
    real(kind=8) r,s
    ! x_coord = (6.5d0+r+5.d0*s)*cos(4.d0*s)
    ! x_coord = (2.d0+r+0.2*sin(5.d0*pi*s))*cos(0.5d0*pi*s)
    x_coord = ((r+1.d0)/2.d0)*cos(3.14159265358979323*(s+1.d0))
  end function x_coord

  real(kind=8) function y_coord(r,s)
    implicit none
    real(kind=8) r,s
    ! y_coord = (6.5d0+r+5.d0*s)*sin(4.d0*s)
    ! y_coord = (2.d0+r+0.2*sin(5.d0*pi*s))*sin(0.5d0*pi*s)
    y_coord = ((r+1.d0)/2.d0)*sin(3.14159265358979323*(s+1.d0))
  end function y_coord
    
end module xycoord
