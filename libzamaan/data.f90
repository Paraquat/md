module dataproc

  use datatypes
  use constants

  implicit none

  contains

! count number of lines in a file
  subroutine countl(filename, n)
    integer                             :: n, ios, i
    character(*)                        :: filename
    character(1)                        :: junk
    integer, parameter                  :: maxrecs = 1000000000

    n = 0
    open(101, file=filename, status='old')
    do i = 1, maxrecs
      read(101,*,end=10) junk
      n = n+1
    end do
    10 close(101)
  end subroutine countl

! Gaussian function
  real(double) function gaussian(a, b, c, x)
    real(double)                        :: a, b, c, x
    gaussian = a*exp(-(x-b)**2/(2*c**2))
  end function gaussian

! Smoothing of histogram via convolution with Gaussian
subroutine smooth_gauss(dataset, sigma, ngrids, gridspacing)
  ! passed variables
  real(double), dimension(:), intent(inout)   :: dataset
  real(double), intent(in)                    :: sigma, gridspacing
  integer, intent(in)                         :: ngrids

  ! local variables
  real(double), dimension(:), allocatable     :: data_temp, temp
  real(double)                                :: e, e1, denom
  integer                                     :: i, j

  allocate(data_temp(ngrids), temp(0:ngrids-1))
  data_temp = zero

  do i=0,ngrids-1
    e = gridspacing*i + gridspacing/two
    temp(i) = gaussian(one, zero, sigma, e)
  end do

  do i=1,ngrids
    e = gridspacing*(i-1) + gridspacing/two
    denom = zero
    do j=1,ngrids
      e1 = gridspacing*(j-1) + gridspacing/two
      data_temp(i) = data_temp(i) + dataset(j) * temp(abs(i-j))
      denom = denom + temp(abs(i-j))
    end do
    dataset(i) = data_temp(i)/denom
  end do

  deallocate(data_temp, temp)

end subroutine smooth_gauss

! integration via trapezium rule
  subroutine trapezium(ydata, xdata, xmin, xmax, binwidth, integral)
    real(double), dimension(:)          :: ydata, xdata
    real(double)                        :: xmin, xmax, binwidth, integral
    integer                             :: nbins, ibin, ixmin, ixmax

    nbins = size(xdata)
    ixmin = int(xmin/binwidth)
    ixmax = int(xmax/binwidth)

    integral = 0.0
    do ibin = ixmin, ixmax-1
      integral = integral + ydata(ibin) + ydata(ibin+1)
    end do
    integral = integral * binwidth / 2.0
  end subroutine trapezium

! cumulative integral via trapezium rule
  subroutine cum_trapezium(ydata, binwidth, cumint)
    real(double), dimension(:)          :: ydata, cumint
    real(double)                        :: binwidth
    integer                             :: nbins, ibin

    nbins = size(ydata)

    cumint = 0.0
    do ibin = 2, nbins-1
      cumint(ibin) = cumint(ibin-1) + ydata(ibin) + ydata(ibin+1)
    end do
    cumint = cumint * binwidth / 2.0
  end subroutine cum_trapezium
end module dataproc
