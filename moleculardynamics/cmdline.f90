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
  character(40)       :: ttype = 'nhc'
  character(40)       :: btype = 'mttk'
  character(40)       :: init_distr = 'uniform'
  character(40)       :: pbc_method = 'mic'
  character(3)        :: ensemble = 'nve'
  real(double)        :: dt
  real(double)        :: T_ext
  real(double)        :: P_ext
  integer             :: nsteps
  integer             :: dump_freq = 1
  integer             :: cp_freq = 1000
  integer             :: n_nhc = 5
  logical             :: restart = .false.
  logical             :: shift = .false.
  logical             :: comv = .false.
  logical             :: cart = .false.
  real(double), dimension(:), allocatable :: nhc_mass
  real(double)        :: boxm = one
  real(double)       :: tau_T = one
  real(double)       :: tau_P = one

  contains
    procedure :: print_help
    procedure :: get_args
end type type_cmdline

contains

subroutine print_help(cmdl)

  class(type_cmdline)           :: cmdl

  write(*,'(a,a,a)') 'Usage: ', trim(cmdl%exe), ' [OPTIONS]'
  write(*,'(a)') 'Will look for md.in if no options are given'
  write(*,*)
  write(*,'(a)') 'Mandatory arguments:'
  write(*,'(a)') '-ns, --nsteps   [nsteps]    number of md steps'
  write(*,'(a)') '-dt, --timestep [dt]        md time step'
  write(*,'(a)') '-T,  --temp     [T_ext]     Initial/target temperature'
  write(*,'(a)')
  write(*,'(a)') 'Optional arguments:'
  write(*,'(a)') '--restart       [restart]   restart from checkpoint'
  write(*,'(a)') '--cp            [cp_freq]   checkpoint frequency (default 1000)'
  write(*,'(a)') '-pp, --ppfile   [ppfile]    input pair potential file (default pp.in)'
  write(*,'(a)') '-i,  --sfile    [sfile]     input structure file (default cell.in)'
  write(*,'(a)') '-o,  --out      [outfile]   output file (default md.out)'
  write(*,'(a)') '-e,  --ensemble [ensemble]  MD ensemble (default nve)'
  write(*,'(a)') '-s,  --shift    [shift]     shift potential/force (default F)'
  write(*,'(a)') '-cv, --comv     [comv]      remove COM velocity (default F)'
  write(*,'(a)') '-c,  --cart     [cart]      Input coordinates are Cartesian (default F)'
  write(*,'(a)') '-d,  --dumpfreq [freq]      Frequency of dump in time steps'
  write(*,'(a)') '-th, --thermo   [ttype]     Thermostat type'
  write(*,'(a)') '-tT, --tau_T    [tau_T]     Frequency of thermostat propagation in time steps'
  write(*,'(a)') '-nh, --n_nhc    [n_nhc]     Length of Nose-Hoover chain (default 5)'
  write(*,'(a)') '-ba, --btype    [btype]     Barostat type'
  write(*,'(a)') '-P,  --press    [P_ext]     External pressure'
  write(*,'(a)') '-Wh, --box_mass [boxm]      Box mass'
  write(*,'(a)') '--pbc           [pbc]       PBC method'

end subroutine print_help

subroutine get_args(cmdl, args)

  class(type_cmdline)           :: cmdl
  logical, intent(out)          :: args

  integer                       :: i, j
  character(40)                 :: arg

  call get_command_argument(0,arg)
  cmdl%exe=arg

  if (command_argument_count() == 0) then
    call cmdl%print_help()
    write(*,"(a)") "No command line arguments."
    args = .false.
  else
    args = .true.
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
      case ('-P', '--press')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%P_ext
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
      case ('-d', '--dumpfreq')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%dump_freq
        i=i+2
      case ('-tT', '--tau_T')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%tau_T
        i=i+2
      case ('-nh', '--n_nhc')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%n_nhc
        i=i+2
        allocate(cmdl%nhc_mass(cmdl%n_nhc))
        cmdl%nhc_mass = one
      case ('-Wh', '--box_mass')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%boxm
        i=i+2
      case ('-th', '--thermo')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%ttype
        i=i+2
      case ('-ba', '--btype')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%btype
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
      case ('--restart')
        cmdl%restart = .true.
        i=i+1
      case ('--cp')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%cp_freq
        i=i+2
      case ('--pbc')
        call get_command_argument(i+1, arg)
        read(arg,*) cmdl%pbc_method
        i=i+2
      case default
        write(*,'(a,a,/)') 'Unrecognised command line option: ', arg
        call cmdl%print_help()
        stop
      end select
    end do
  end if

end subroutine get_args

end module cmdline
