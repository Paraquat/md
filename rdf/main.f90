program rdf

use datatypes
use progress
use pairdist
use dataproc
use cmdline

implicit none

type(type_pairdist)                 :: pd
type(type_cmdline)                  :: opts
type(type_cell)                     :: p

character(40)                       :: outfile, junk
integer                             :: i, j, k, cellrad, nframes, ios, iounit

real(double)                        :: cutoff
real(double), dimension(:), allocatable :: b    ! coherent neutron scattering length

! for nearest neigbour of atom with label centreid
integer                             :: nneigh

call opts%get_args()

iounit = 101
cellrad = 1
b = (/1/)

open(iounit, file=opts%infile, status='old', iostat=ios)
read(iounit,*) junk, nframes
read(iounit,*) junk
read(iounit,*) junk
read(iounit,*) p%h(1,:)
read(iounit,*) p%h(2,:)
read(iounit,*) p%h(3,:)
read(iounit,*) junk
read(iounit,*) p%nat, k

allocate(p%species(p%nat), p%r(p%nat,3), p%rcart(p%nat,3))
call p%invert_lat()

do j=1,p%nat
  read(iounit,*) p%species(j), p%rcart(j,:)
  p%r(j,:) = p%cart2frac(p%rcart(j,:))
end do
call p%count_species()
call p%cut_ortho(cutoff)
call pd%init_pd(p, cellrad, cutoff, opts%gwidth, opts%binwidth, opts%rmin)
pd%smooth_on = opts%smooth
allocate(pd%p%dt(pd%p%nat,pd%p%nat))

if (opts%noH .eqv. .true.) then
  pd%ignorespec = 'H'
end if
call pd%p%get_dt()
call pd%update_rdist_dt()

do i=2,nframes
  read(iounit,*) junk
  read(iounit,*) junk
  read(iounit,*) pd%p%h(1,:)
  read(iounit,*) pd%p%h(2,:)
  read(iounit,*) pd%p%h(3,:)
  read(iounit,*) junk
  read(iounit,*) junk

  call pd%p%invert_lat()

  do j=1,pd%p%nat
    read(iounit,*) pd%p%species(j), pd%p%rcart(j,:)
    pd%p%r(j,:) = pd%p%cart2frac(pd%p%rcart(j,:))
  end do
  call p%cut_ortho(cutoff)

  call pd%p%get_dt()
  call pd%update_rdist_dt()
end do

if (opts%freq .eqv. .true.) then
  call pd%get_cumfreq()
  call pd%write_count(opts%outfile, opts%outfile_int)
else
  call pd%norm_rdist(nframes)
  call pd%write_gr(opts%outfile)
end if
close(iounit)

end program rdf
