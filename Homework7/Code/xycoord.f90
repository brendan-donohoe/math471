module xycoord
  real(kind = 8), parameter :: pi = acos(-1.d0)
contains
  
  real(kind=8) function x_coord(r,s)
    implicit none
    real(kind=8) r,s
    x_coord = (2.d0+r+0.4d0*sin(s))*cos(s) 
    !x_coord = r
    !x_coord = r+0.1d0*s
  end function x_coord

  real(kind=8) function y_coord(r,s)
    implicit none
    real(kind=8) r,s
    y_coord = (2.d0+r+0.4d0*sin(s))*sin(s) 
    !y_coord = s
    !y_coord = s + s*r**2
  end function y_coord
    
end module xycoord
