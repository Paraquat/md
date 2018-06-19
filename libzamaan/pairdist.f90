module pairdist

use constants
use progress
use cell
use dataproc
use fft
use datatypes

implicit none

public :: type_pairdist

type type_pairdist

  type(type_cell)                         :: p

  integer                                 :: maxspecies=5
  integer                                 :: nat, nbins, nsteps, nqpt, nframes
  character(2)                            :: ignorespec

  real(double)                            :: gwidth, rho, vol, cut, cutsq, &
                                             delr, delq, sqconst, rmin, &
                                             qmin, qmax
  real(double), dimension(3,3)            :: h

  real(double), allocatable, dimension(:) :: bins, nfac_total, &
                                             gr_total, sq_total, sq_sum
  real(double), allocatable, dimension(:,:,:) :: nfac, gr, sq
  integer, allocatable, dimension(:)      :: total, cum_total
  integer, allocatable, dimension(:,:,:)  :: freq, cum_freq
  character(2), allocatable, dimension(:) :: sname
  integer, allocatable, dimension(:)      :: scount
  logical                                 :: smooth_on = .true.

  contains
    procedure :: init_pd
    procedure :: reset_pd
    procedure :: update_rdist
    procedure :: get_neighbours
    procedure :: get_nfac
    procedure :: smooth_gr
    procedure :: norm_rdist
    procedure :: get_cumfreq
    procedure :: get_sfac_grdft
    procedure :: get_sfac_grfft
    procedure :: sum_sfac
    procedure :: write_gr
    procedure :: write_sq
    procedure :: write_count
end type type_pairdist

contains

! initialise the pair distribution object from a cell object
! gwidth = width for gaussian smoothing
! delr = width of histogram bins for r
subroutine init_pd(pd, cell, cut, gwidth, delr, rmin)

  ! passed variables
  class(type_pairdist), intent(inout)     :: pd
  type(type_cell), intent(in)             :: cell
  type(type_cell)                         :: p
  real(double), intent(in)                :: cut, gwidth, delr, rmin

  ! local variables
  real(double)                            :: min_side

  pd%p = cell
  pd%h = p%h
  pd%cut = cut                        ! distance cutoff for g(r)
  pd%cutsq = cut**2                   ! distance cutoff squared
  pd%gwidth = gwidth
  pd%delr = delr                      ! bin width in real space
  pd%nbins = ceiling(pd%cut/pd%delr)  ! number of real space bins)
  pd%nsteps = 0                       ! number of MD steps
  pd%qmin = 2*pi/cut
  pd%rmin = rmin                      ! minimum value of r
  pd%ignorespec = ''                  ! species to omit from g_total(r)
  pd%nframes = 0                      ! number for frames for normalisation

  min_side = min(cell%h(1,1), cell%h(2,2), cell%h(3,3))
  if (min_side > pd%cut) then
    write(*,'(2x,a)') "Warning: rdf cutoff is less than shortest box side"
    write(*,'(2x,a,f8.4)') "Adjusting rdfcut to ", min_side
    pd%cut = min_side
  end if

  ! constants
  pd%vol = pd%p%volume()                      ! volume
  pd%rho = real(pd%p%nat,double)/pd%vol              ! density

  allocate(pd%bins(pd%nbins))
  allocate(pd%total(pd%nbins))
  allocate(pd%cum_total(pd%nbins))
  allocate(pd%nfac_total(pd%nbins))
  allocate(pd%gr_total(pd%nbins))
  allocate(pd%freq(pd%nbins,pd%p%nspec,pd%p%nspec))
  allocate(pd%cum_freq(pd%nbins,pd%p%nspec,pd%p%nspec))
  allocate(pd%nfac(pd%nbins,pd%p%nspec,pd%p%nspec))
  allocate(pd%gr(pd%nbins,pd%p%nspec,pd%p%nspec))

  pd%total = 0
  pd%cum_total = 0
  pd%freq = 0
  pd%cum_freq = 0

end subroutine init_pd

subroutine reset_pd(pd)

  class(type_pairdist), intent(inout)     :: pd

  pd%total = 0
  pd%cum_total = 0
  pd%freq = 0
  pd%cum_freq = 0
end subroutine reset_pd

! Update pair distribution using distance table
subroutine update_rdist(pd, p)

  class(type_pairdist), intent(inout)     :: pd
  type(type_cell), intent(in)             :: p

  integer                                 :: iat, jat, ispec, jspec, ig

  ! bin atomic separations
  do iat=1,p%nat
    do jat=iat+1,p%nat
      if (p%dt(iat,jat) .lt. pd%cut) then
        if (p%species(iat) .ne. pd%ignorespec .and. p%species(jat) .ne. pd%ignorespec) then
          ig = int((p%dt(iat,jat)+pd%delr)/pd%delr)
          pd%total(ig) = pd%total(ig)+2
          do ispec=1,p%nspec
            do jspec=ispec,p%nspec
              if (p%species(iat) .eq. p%spec(ispec) .and. &
                  p%species(modulo(jat,p%nat)+1) .eq. p%spec(jspec)) then
                pd%freq(ig,ispec,jspec) = pd%freq(ig,ispec,jspec)+2
              end if
            end do
          end do
        end if
      end if
    end do
  end do
  pd%nframes = pd%nframes + 1

end subroutine update_rdist

! Use distance table to get list of neighbours of atom centreid with a
! given radius
subroutine get_neighbours(pd, centreid, radius, nneigh, nlist)

  class(type_pairdist), intent(inout)     :: pd

  integer, intent(in)                     :: centreid
  real(double), intent(in)                :: radius
  integer, intent(out)                    :: nneigh
  integer, dimension(:), intent(out)      :: nlist

  integer                                 :: i

  nneigh=0
  do i=1,pd%p%nat
    if (pd%p%dt(centreid,i) .lt. radius) then
      nneigh=nneigh+1
      nlist(nneigh)=i
    end if
  end do

end subroutine get_neighbours

! Compute the ideal gas normalisation factors
subroutine get_nfac(pd, nframes)

  class(type_pairdist), intent(inout)     :: pd

  integer, intent(in)                     :: nframes

  real(double)                            :: vshell, grconst1, grconst2, &
                                             rho_reduced
  integer                                 :: i, ispec, jspec, nat_reduced

  rho_reduced = pd%rho
  nat_reduced = pd%p%nat
  if (pd%ignorespec .ne. '') then
    do i=1,pd%p%nspec
      if (pd%p%spec(i) .eq. pd%ignorespec) then
        ispec = i
        exit
      end if
    end do
    nat_reduced = pd%p%nat-pd%p%spec_count(ispec)
    rho_reduced = real(nat_reduced,double)/pd%vol
  end if
  grconst1 = rho_reduced*nat_reduced
  do i=1,pd%nbins
    vshell = 4.*pi*((i+1)**3 - i**3)*pd%delr**3/3.
    pd%nfac_total(i) = vshell*grconst1*nframes
    do ispec=1,pd%p%nspec
      do jspec=ispec,pd%p%nspec
        grconst2 = pd%rho*pd%p%spec_count(ispec)*pd%p%spec_count(jspec)/pd%p%nat
        pd%nfac(i,ispec,jspec) = vshell*grconst2*nframes
      end do
    end do
  end do
end subroutine get_nfac

! Smooth g(r)
subroutine smooth_gr(pd)

  class(type_pairdist), intent(inout)     :: pd

  integer, parameter                      :: lwindow = 10
  integer                                 :: i,j

  call smooth_gauss(pd%gr_total,pd%gwidth,pd%nbins,pd%delr)
  do i=1,pd%p%nspec
    do j=i,pd%p%nspec
      call smooth_gauss(pd%gr(:,i,j),pd%gwidth,pd%nbins,pd%delr)
    end do
  end do
end subroutine smooth_gr

! Normalise the pair distribution to get g(r)
subroutine norm_rdist(pd)

  class(type_pairdist), intent(inout)     :: pd

  integer, parameter                      :: lwindow=10
  integer                                 :: i, j

  call pd%get_nfac(pd%nframes)

  pd%gr_total = pd%total/pd%nfac_total
  pd%gr = pd%freq/pd%nfac

  if (pd%smooth_on .eqv. .true.) then
    call pd%smooth_gr()
  end if

end subroutine norm_rdist

! Compute cumulative frequency of distances
subroutine get_cumfreq(pd)

  class(type_pairdist), intent(inout)     :: pd

  integer                                 :: i

  pd%cum_total(1) = pd%total(1)
  pd%cum_freq(1,:,:) = pd%freq(1,:,:)

  do i=2,pd%nbins
    pd%cum_total(i) = pd%cum_total(i-1) + pd%total(i)
    pd%cum_freq(i,:,:) = pd%cum_freq(i-1,:,:) + pd%freq(i,:,:)
  end do

end subroutine get_cumfreq

! Calculate the structure factor via slow discrete Fourier transform of g(r)
subroutine get_sfac_grdft(pd,nqpt,qmax)

  class(type_pairdist), intent(inout)     :: pd

  integer                                 :: ispec, jspec, i, nqpt
  real(double)                            :: delq, qmax

  pd%qmax=qmax
  pd%nqpt=nqpt
  if (allocated(pd%sq_total) .eqv. .false.) then
    allocate(pd%sq_total(nqpt))
    allocate(pd%sq(pd%nqpt,pd%p%nspec,pd%p%nspec))
  end if

  call grdft(pd%gr_total,pd%sq_total,qmax,pd%delr,pd%delq)
  do ispec=1,pd%p%nspec
    do jspec=ispec,pd%p%nspec
      call grdft(pd%gr(:,ispec,jspec),pd%sq(:,ispec,jspec),qmax,pd%delr, &
                 pd%delq)
    end do
  end do

end subroutine get_sfac_grdft

! Calculate the structure factor via Fourier transform of g(r)
! Use grfft subroutine in grfft.f90
! THIS IS BROKEN
subroutine get_sfac_grfft(pd,nqpt)

  class(type_pairdist), intent(inout)     :: pd

  integer                                 :: ispec, jspec, i, nqpt
  real(double)                            :: delq

  pd%nqpt=nqpt
  allocate(pd%sq_total(pd%nqpt))
  allocate(pd%sq(pd%nqpt,pd%p%nspec,pd%p%nspec))
  pd%sqconst = 4.*pi*pd%rho*pd%delr**2    ! structure factor constant

  call grfft(pd%gr_total,pd%sq_total,pd%nbins,pd%nqpt,pd%sqconst)
end subroutine get_sfac_grfft

! Compute the total structure factor from partial structure factors
subroutine sum_sfac(pd, b)

  class(type_pairdist), intent(inout)     :: pd

  real(double), dimension(:), intent(in)  :: b

  integer                                 :: iat,ispec,jspec
  real(double)                            :: c1,c2

  if (allocated(pd%sq_sum) .eqv. .false.) then
    allocate(pd%sq_sum(pd%nqpt))
  end if
  pd%sq_sum = 0.

  if (pd%p%nspec .gt. 1) then
    do iat=1,pd%p%nat
      do ispec=1,pd%p%nspec
        do jspec=ispec,pd%p%nspec
          c1 = pd%p%spec_count(ispec)/pd%p%nat
          c2 = pd%p%spec_count(jspec)/pd%p%nat
          pd%sq_sum(iat) = pd%sq_sum(iat) + &
                          b(ispec)*b(jspec)*c1*c2*pd%sq(iat,ispec,jspec)
        end do
      end do
    end do
  else
    pd%sq_sum = pd%sq(:,1,1)
  end if

end subroutine sum_sfac

! Write g(r) to filename
subroutine write_gr(pd, filename)

  class(type_pairdist), intent(inout)     :: pd

  character(40), intent(in)               :: filename

  character(4)                            :: pair
  integer                                 :: fh, iat, ispec, jspec
  real(double)                            :: r

  fh = 104
  open(fh, file=filename)

  write(fh,'(a12)',advance="no") "r"
  write(fh,'(a12)',advance="no") "total"
  do ispec=1,pd%p%nspec
    do jspec=ispec,pd%p%nspec
      pair = trim(pd%p%spec(ispec)) // trim(pd%p%spec(jspec))
      if (ispec .eq. pd%p%nspec .and. jspec .eq. pd%p%nspec) then
        write(fh,'(a12)') pair
      else
        write(fh,'(a12)',advance="no") pair
      end if
    end do
  end do

  do iat=1,pd%nbins
!    write(fh,'(f12.6)',advance="no") pd%delr*(iat-0.5)
    r=pd%delr*iat
    if (r+0.000001 .gt. pd%rmin) then
      write(fh,'(f12.6)',advance="no") r
      write(fh,'(f12.6)',advance="no") pd%gr_total(iat)
      do ispec=1,pd%p%nspec
        do jspec=ispec,pd%p%nspec
          if (ispec .eq. pd%p%nspec .and. jspec .eq. pd%p%nspec) then
            write(fh,'(f12.6)') pd%gr(iat,ispec,jspec)
          else
            write(fh,'(f12.6)',advance="no") pd%gr(iat,ispec,jspec)
          end if
        end do
      end do
    end if
  end do

  close(fh)
end subroutine write_gr

! Write s(q) to filename
subroutine write_sq(pd, filename)

  class(type_pairdist), intent(inout)     :: pd

  character(40), intent(in)               :: filename

  character(4)                            :: pair
  integer                                 :: fh, iq, ispec, jspec
  real(double)                            :: q

  fh = 105
  open(fh, file=filename)

  write(fh,'(a12)',advance="no") "q"
  write(fh,'(a12)',advance="no") "total"
  do ispec=1,pd%p%nspec
    do jspec=ispec,pd%p%nspec
      pair = trim(pd%p%spec(ispec)) // trim(pd%p%spec(jspec))
      if (ispec .eq. pd%p%nspec .and. jspec .eq. pd%p%nspec) then
        write(fh,'(a12)') pair
      else
        write(fh,'(a12)',advance="no") pair
      end if
    end do
  end do

  do iq=1,pd%nqpt
    q = iq*pd%delq
    if (q .gt. pd%qmin .and. q .le. pd%qmax) then
      write(fh,'(f12.6)',advance="no") q
      write(fh,'(f12.6)',advance="no") pd%sq_total(iq)**2
      if (allocated(pd%sq_sum) .eqv. .true.) then
        write(fh,'(f12.6)',advance="no") pd%sq_sum(iq)**2
      end if
      do ispec=1,pd%p%nspec
        do jspec=ispec,pd%p%nspec
          if (ispec .eq. pd%p%nspec .and. jspec .eq. pd%p%nspec) then
            write(fh,'(f12.6)') pd%sq(iq,ispec,jspec)
          else
            write(fh,'(f12.6)',advance="no") pd%sq(iq,ispec,jspec)**2
          end if
        end do
      end do
    end if
  end do

  close(fh)
end subroutine write_sq

! Write distance frequencies to filename
subroutine write_count(pd, filename, filename_int)

  class(type_pairdist), intent(inout)     :: pd

  character(40), intent(in)               :: filename, filename_int

  character(4)                            :: pair
  integer                                 :: fh1, fh2, iat, ispec, jspec

  call pd%get_cumfreq()

  fh1 = 104
  fh2 = 105
  open(fh1, file=filename)
  open(fh2, file=filename_int)

  write(fh1,'(a10)',advance="no") "r"
  write(fh1,'(a12)',advance="no") "total"
  write(fh2,'(a10)',advance="no") "r"
  write(fh2,'(a12)',advance="no") "total"
  do ispec=1,pd%p%nspec
    do jspec=ispec,pd%p%nspec
      pair = trim(pd%p%spec(ispec)) // trim(pd%p%spec(jspec))
      if (ispec .eq. pd%p%nspec .and. jspec .eq. pd%p%nspec) then
        write(fh1,'(a12)') pair
        write(fh2,'(a12)') pair
      else
        write(fh1,'(a12)',advance="no") pair
        write(fh2,'(a12)',advance="no") pair
      end if
    end do
  end do

  do iat=1,pd%nbins
    write(fh1,'(f10.6)',advance="no") pd%delr*(iat-0.5)
    write(fh1,'(i12)',advance="no") pd%total(iat)
    write(fh2,'(f10.6)',advance="no") pd%delr*(iat-0.5)
    write(fh2,'(i12)',advance="no") pd%cum_total(iat)
    do ispec=1,pd%p%nspec
      do jspec=ispec,pd%p%nspec
        if (ispec .eq. pd%p%nspec .and. jspec .eq. pd%p%nspec) then
          write(fh1,'(i12)') pd%freq(iat,ispec,jspec)
          write(fh2,'(i12)') pd%cum_freq(iat,ispec,jspec)
        else
          write(fh1,'(i12)',advance="no") pd%freq(iat,ispec,jspec)
          write(fh2,'(i12)',advance="no") pd%cum_freq(iat,ispec,jspec)
        end if
      end do
    end do
  end do

  close(fh1)
end subroutine write_count

end module pairdist
