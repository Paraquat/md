module rng

! Generate random numbers using intrinsic compiler function. Sets seed using
! wall time.

use datatypes

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

end module rng
