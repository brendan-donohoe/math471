module xycoord
  real(kind = 8), parameter :: pi = acos(-1.d0)
  save
contains
  
  real(kind=8) function x_coord(r,s)
    implicit none
    real(kind=8) r,s
    
    ! x_coord = (2.d0+r+0.2*sin(5.d0*pi*s))*cos(0.5d0*pi*s)

    ! FUNCTION 0:
    x_coord = ((r+1.d0)/2.d0)*cos(pi*(s+1.d0))

    ! FUNCTION 1:
    ! x_coord = r + 0.1*s

    ! FUNCTION 2:
    ! x_coord = (2.d0+r+0.1*sin(5.d0*pi*s))*cos(0.5d0*pi*s)

    ! FUNCTION 3:
    ! x_coord = r
  end function x_coord

  real(kind=8) function y_coord(r,s)
    implicit none
    real(kind=8) r,s

    ! y_coord = (2.d0+r+0.2*sin(5.d0*pi*s))*sin(0.5d0*pi*s)

    ! FUNCTION 0:
    y_coord = ((r+1.d0)/2.d0)*sin(pi*(s+1.d0))

    ! FUNCTION 1:
    ! y_coord = s

    ! FUNCTION 2:
    ! y_coord = (2.d0+r+0.1*sin(5.d0*pi*s))*sin(0.5d0*pi*s)

    ! FUNCTION 3:
    ! y_coord = s + s*r**2
  end function y_coord
    
end module xycoord
