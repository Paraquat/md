module vector

  use constants
  use datatypes

  implicit none

  contains

  ! Calculate length of a three vector
  function modulus(vector) result(l)
    real(double), dimension(3)  :: vector
    real(double)                :: l

    l = sqrt(sum(vector**2))
  end function modulus

  ! Normalise a vector
  function norm(vector) result(n)
    real(double), dimension(3)  :: vector, n
    real(double)                :: l

    l = modulus(vector)
    n = vector / l
  end function norm

  ! Cross product
  subroutine cross_product(u, v, cross)
    real(double), dimension(3)  :: u, v, cross

    cross(1) = u(2)*v(3) - u(3)*v(2)
    cross(2) = u(3)*v(1) - u(1)*v(3)
    cross(3) = u(1)*v(2) - u(2)*v(1)
  end subroutine cross_product

  ! angle between vectors u and v in degrees
  subroutine vangle(u, v, thetad)
    real(double), dimension(3)  :: u, v
    real(double)                :: angle, modu, modv, thetad, dp, prod

    modu = sqrt(sum(u**2))
    modv = sqrt(sum(v**2))
    dp = sum(u*v)
    prod = dp/(modu*modv)
    if (abs(prod) .gt. 1.0) then
      if ((abs(prod) - 1.0) .lt. zero) then
        prod = -1.0
      else
        write(*,*) "Trying to calculate acos(x) where abs(x) > 1"
        stop
      end if
    end if
    thetad = rad2deg*acos(prod)
  end subroutine vangle

end module vector
