program rdf

use progress
use pairdist
use dataproc
use cmdline

implicit none

type(type_pairdist)                 :: pd
type(type_cmdline)                  :: opts
type(type_cell)                     :: p

character(40)                       :: infile, outfile
integer                             :: i, j, nspec, cellrad, nframes

real                                :: cutoff, gwidth
real, dimension(:), allocatable     :: b    ! coherent neutron scattering length

! for nearest neigbour of atom with label centreid
integer                             :: nneigh

call opts%get_args()

cellrad = 1
gwidth = 0.1
nframes = 1
b = (/1/)

if (opts%lammps .eqv. .true.) then
  call p%read_lammps(opts%infile, opts%nspec, opts%speclist)
!else if (opts%outcar == .true.) then
!  pass
else if (opts%xsf .eqv. .true.) then
  call p%read_xsf(opts%infile)
end if
call p%cut_ortho(cutoff)
call pd%init_pd(p, cellrad, cutoff, gwidth, opts%binwidth, opts%rmin)

pd%smooth_on = opts%smooth
if (opts%noH .eqv. .true.) then
  pd%ignorespec = 'H'
end if
call pd%get_dt()
!call pd%update_rdist()
call pd%update_rdist_dt()
if (opts%freq .eqv. .true.) then
  call pd%get_cumfreq()
  call pd%write_count(opts%outfile, opts%outfile_int)
else
  call pd%norm_rdist(nframes)
  call pd%write_gr(opts%outfile)
end if

end program rdf
