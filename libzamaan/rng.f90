module rng

! Generate random numbers using intrinsic compiler function. Sets seed using
! wall time.

use datatypes
use constants

implicit none

contains

! seed the generator using wall time
subroutine init_rand()
  integer, parameter                  :: maxr = huge(0) - 1
  real(double)                        :: rr
  integer                             :: i, j, n, pid, c, cr, cmax, iters
  integer, dimension(8)               :: time
  integer, dimension(3)               :: v
  integer, allocatable                :: seed(:)

  pid = getpid()
  call random_seed(size=n)
  allocate(seed(n))
  call random_seed()

  do j = 1, n
    call system_clock(c, cr, cmax)
    iters = int(mod(c, pid), 4)
    do i = 1, iters + 1
      call random_number(rr)
    end do
    seed(j) = int(rr * maxr)
  end do

  call random_seed(put=seed)
  deallocate(seed)
end subroutine init_rand

! random real number between 0 and 1
subroutine rand(rr)
  real(double)                        :: rr

  call random_number(rr)
end subroutine rand

! random integer ri between rmin and rmax
subroutine randint(rmin, rmax, ri)
  integer                             :: ri, rmin, rmax
  real(double)                        :: rr

  call rand(rr)
  ri = floor(real(rmax + 1 - rmin, double) * rr) + rmin
end subroutine randint

! two random number generated from a Gaussian distribution
subroutine boxmuller(sigma, mu, y1, y2)
  real(double), intent(in)  :: sigma, mu
  real(double), intent(out) :: y1, y2
  real(double)              :: x1, x2

  call rand(x1)
  call rand(x2)
  y1 = sqrt(-two*log(x1))*cos(twopi*x2)
  y2 = sqrt(-two*log(x1))*sin(twopi*x2)
end subroutine boxmuller

! polar version of Box-Muller (no sine/cosine calls)
subroutine boxmuller_polar(sigma, mu, y1, y2)
  real(double), intent(in)  :: sigma, mu
  real(double), intent(out) :: y1, y2
  real(double)              :: x1, x2, w

  do
    call rand(x1)
    call rand(x2)
    y1 = two*x1-one
    y2 = two*x2-one
    w = y1**2 + y2**2
    if (w < 1.0) exit
  end do
  w = sqrt(-two*log(w)/w)
  y1 = y1*w*sigma + mu
  y2 = y2*w*sigma + mu

end subroutine boxmuller_polar

end module rng
