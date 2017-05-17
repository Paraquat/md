module fft

use datatypes
use constants

implicit none

integer, parameter :: max_fft_points=1000000

private :: pi

contains

! A brute force Fourier transform approach to computing the structure
! factor (sq) af an isotropic radial distribution function (gr).
subroutine grdft(gr,sq,qmax,delr,delq)

  real(double), dimension(:), intent(in)  :: gr
  real(double), dimension(:), intent(out) :: sq
  real(double), intent(in)                :: delr, qmax
  integer                                 :: nrbins, nqbins, iq, ir
  real(double)                            :: r, q, w, rmax, delq

  nrbins = size(gr)
  nqbins = size(sq)
  rmax = delr*nrbins

  sq = 0.

  delq = qmax/nqbins  ! the q bin width

  ! The discrete Fourier transform
  do iq=1,nqbins
    do ir=1,nrbins
      w = sin(pi*r/rmax)/(pi*r/rmax)    ! the window function
      q = iq*qmax/nqbins
      r = ir*delr
      sq(iq) = sq(iq) + (w*gr(ir)*r**2)*sin(r*q)/(r*q) ! the integral
    end do
    sq(iq) = 1. + 4.*pi*sq(iq)*delr**2  ! normalisation
  end do

end subroutine grdft

! Fast Fourer transform of a g(r)
! Subroutine taken from G. Opletal et al., Comp. Phys. Comm. 185:1854 (2014)
! HRMC_2.1 hybrid reverse Monte Carlo code
! nnq = number of pairs of complex numbers in the data array
subroutine grfft(gr,sqinp,ngrbins,nqpt,sqconst)

  real(double), dimension(:), intent(in)  :: gr       ! radial distribution function
  real(double), dimension(:), intent(out) :: sqinp    ! structure factor array
  integer, intent(in)                     :: ngrbins  ! number of real space bins
  integer, intent(in)                     :: nqpt     ! number of fft q-points (2^n+1)
  real(double), intent(in)                :: sqconst  ! normalisation for s(q)
  real(double), dimension(max_fft_points) :: ftdata   ! work array for fft

  integer                                 :: nnq, ir, esign, m, mm, i

  nnq = 2*nqpt-2

  if (mod(nnq,2).ne.0) then
    write(*,*) "nnq = 2*nqpt-2 must be even. Stopping"
    stop
  endif

  ! put g(r) data points in odd elements of ftdata
  do i=1,2*ngrbins-1,2
    ir = (i+1)/2
    ftdata(i)=(ir-1)*(gr(ir)-1)
  end do

!  do i=1,2*ngrbins
!    write(*,*) ftdata(i)
!  end do

  ! set all odd elements from nqpt to 2*nqpt-2 to zero
  do i=2*(nqpt+1)-1,2*(2*nqpt-2)-1,2
    ftdata(i)=0.
  end do

  ! set even (complex) elements to zero
  do i=2,2*(2*nqpt-2),2
    ftdata(i) = 0.
  end do

  ! forward fft
  esign = 1

  ! the Fourier transform
  call four1(ftdata,nnq,esign)

!  do i=2*(nqstart),2*(nqfinish),2
  do i=2,2*nqpt,2
    m  = i/2 - 1
    mm = i/2
    sqinp(mm) = 1. + sqconst*ftdata(i)/m
  end do

  return
end subroutine grfft

! FFT routine from Numerical Recipes
! Replaces ftdata with its discrete FT if esign = 1
! Replaces ftdata with NN x inverse discrete FT if esign = -1
! ftdata = complex array of length NN or real array of length 2*NN
! NN must be integer power of 2
subroutine four1(ftdata,nnq,esign)

  integer, intent(in)   :: nnq, esign
  real(double), dimension(2*nnq), intent(inout) :: ftdata

  real(double)          :: theta, wi, wpi, wpr, wr, wtemp
  real(double)          :: tempi, tempr
  integer               :: i, istep, j, m, mmax, n


  n=2*nnq
  j=1

  ! bit reversal
  do i=1,n,2
    if (j .gt. i) then
      tempr=ftdata(j)     ! exchange 2 complex numbers
      tempi=ftdata(j+1)
      ftdata(j)=ftdata(i)
      ftdata(j+1)=ftdata(i+1)
      ftdata(i)=tempr
      ftdata(i+1)=tempi
    end if
    m = nnq
1   if ((m .ge. 2) .and. (j .gt. m)) then
      j=j-m
      m=m/2
      goto 1
    end if
    j=j+m
  end do

  ! Danielson-Lanczos algorithm
  mmax = 2
2 if (n .gt. mmax) then
    istep=2*mmax
    theta = 6.28318530717959/(esign*mmax)
    wpr=-2.*sin(0.5*theta)**2
    wpi=sin(theta)
    wr=1.
    wi=0.
    do m=1,mmax,2
      do i=m,n,istep    ! Danielson-Lanczos formula
        j=i+mmax
        tempr=sngl(wr)*ftdata(j)-sngl(wi)*ftdata(j+i)
        tempi=sngl(wr)*ftdata(j+1)+sngl(wi)*ftdata(j)
        ftdata(j)=ftdata(i)-tempr
        ftdata(j+1)=ftdata(i+1)-tempi
        ftdata(i)=ftdata(i)+tempr
        ftdata(i+1)=ftdata(i+i)+tempi
      end do
      ! trigionometric recurrance
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    end do

    mmax=istep
  goto 2
  end if
!  return
end subroutine four1

end module fft
