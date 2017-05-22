module cmdline

use datatypes
use constants

implicit none

type type_cmdline

  character(40)       :: exe
  integer             :: nargs

  character(40)       :: ppfile = 'pp.in'
  character(40)       :: sfile = 'cell.in'
  character(40)       :: outfile = 'md.out'
  character(3)        :: ensemble = 'nve'
  real(double)        :: dt
  real(double)        :: T_ext
  integer             :: nsteps
  logical             :: shift = .false.
  logical             :: comv = .false.
  logical             :: cart = .false.

  contains
    procedure :: print_help
    procedure :: get_args
end type type_cmdline

contains

subroutine print_help(cmdl)

  class(type_cmdline)           :: cmdl

  write(*,'(a,a,a)') 'Usage: ', trim(cmdl%exe), ' [OPTIONS]'
  write(*,*)
  write(*,'(a)') 'Mandatory arguments:'
  write(*,'(a)') '-ns, --nsteps   [nsteps]    number of md steps'
  write(*,'(a)') '-dt, --timestep [dt]        md time step'
  write(*,'(a)') '-T,  --temp     [T_ext]     Initial/target temperature'
  write(*,'(a)')
  write(*,'(a)') 'Optional arguments:'
  write(*,'(a)') '-pp, --ppfile   [ppfile]    input pair potential file (default pp.in)'
  write(*,'(a)') '-i,  --sfile    [sfile]     input structure file (default cell.in)'
  write(*,'(a)') '-o,  --out      [outfile]   output file (default md.out)'
  write(*,'(a)') '-e,  --ensemble [ensemble]  MD ensemble (default nve)'
  write(*,'(a)') '-s,  --shift    [shift]     shift potential/force (default F)'
  write(*,'(a)') '-rv, --comv     [comv]      remove COM velocity (default F)'
  write(*,'(a)') '-c,  --cart     [cart]      Input coordinates are Cartesian (default F)'

end subroutine print_help

subroutine get_args(cmdl)

  class(type_cmdline)           :: cmdl

  integer                       :: i, j
  character(40)                 :: arg

  call get_command_argument(0,arg)
  cmdl%exe=arg

  if (command_argument_count() == 0) then
    call cmdl%print_help()
    stop "No command line arguments."
  end if

  i=1
  do while (i .le. command_argument_count())
    call get_command_argument(i, arg)

    select case (arg)
    case ('-h', '--help')
      call cmdl%print_help()
      stop
    case ('-ns', '--nsteps')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%nsteps
      i=i+2
    case ('-dt', '--timestep')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%dt
      i=i+2
    case ('-T', '--temp')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%T_ext
      i=i+2
    case ('-pp', '--ppfile')
      call get_command_argument(i+1, arg)
      cmdl%ppfile=trim(arg)
      i=i+2
    case ('-i', '--sfile')
      call get_command_argument(i+1, arg)
      cmdl%sfile=trim(arg)
      i=i+2
    case ('-o', '--out')
      call get_command_argument(i+1, arg)
      cmdl%outfile=trim(arg)
      i=i+2
    case ('-e', '--ensemble')
      call get_command_argument(i+1, arg)
      cmdl%ensemble=trim(arg)
      i=i+2
    case ('-s', '--shift')
      cmdl%shift=.true.
      i=i+1
    case ('-cv', '--comv')
      cmdl%comv=.true.
      i=i+1
    case ('-c', '--cart')
      cmdl%cart=.true.
      i=i+1
    case default
      write(*,'(a,a,/)') 'Unrecognised command line option: ', arg
      call cmdl%print_help()
      stop
    end select
  end do

end subroutine get_args

end module cmdline
