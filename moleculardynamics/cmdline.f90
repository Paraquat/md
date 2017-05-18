module cmdline

implicit none

type type_cmdline

  character(40)       :: exe
  integer             :: nargs

  character(40)       :: infile = 'md.in'
  character(40)       :: outfile = 'md.out'
  real                :: cutoff
  real                :: binwidth = 0.01
  real                :: gwidth = 0.1
  real                :: rmin = 0.0
  integer             :: nspec
  logical             :: smooth = .false.
  logical             :: freq = .false.
  logical             :: noH = .false.
  logical             :: xsf = .false.
  logical             :: outcar = .false.
  logical             :: lammps = .false.
  character(2), dimension(10) :: speclist = 'a'

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
  write(*,'(a)') '-ns, --nspec    [nspec]     number of atomic species'
  write(*,'(a)') '-s,  --species  [C H O...]  atomic species ordered by label'
  write(*,'(a)')
  write(*,'(a)') 'Optional arguments:'
  write(*,'(a)') '-i,  --inp      [infile]    input structure file (lammps format)'
  write(*,'(a)') '-o,  --out      [outfile]   output g(r) file'
  write(*,'(a)') '-o2,  --out2    [outfile2]  output integrated g(r) file'
  write(*,'(a)') '-c,  --cut      [cutoff]    distance cutoff for g(r)'
  write(*,'(a)') '-dr, --binwidth [binwidth]  histogram bin width for g(r)'
  write(*,'(a)') '-rm, --rmin     [rmin]      minimum value of r'
  write(*,'(a)') '-sm, --smooth               switches Gaussian smoothing on'
  write(*,'(a)') '-f,  --freq     [freq]      compute frequencies instead of g(r)'
  write(*,'(a)') '-H, --noH                   omit H from total g(r)'
  write(*,'(a)') '-g, --gwidth   [gwidth]     window size for Gaussian smoothing'
  write(*,'(a)') '--xsf          [xsf]        read a .xsf trajectory'
  write(*,'(a)') '--outcar       [outcar]     read a vasp trajectory'
  write(*,'(a)') '--lammps       [lammps]     read a lammps trajectory'

end subroutine print_help

subroutine get_args(cmdl)

  class(type_cmdline)           :: cmdl

  integer                       :: i, j
  character(40)                 :: arg

  call get_command_argument(0,arg)
  cmdl%exe=arg

  i=1
  do while (i .le. command_argument_count())
    call get_command_argument(i, arg)

    select case (arg)
    case ('-h', '--help')
      call cmdl%print_help()
      stop
    case ('-i', '--inp')
      call get_command_argument(i+1, arg)
      cmdl%infile=trim(arg)
      i=i+2
    case ('-o', '--out')
      call get_command_argument(i+1, arg)
      cmdl%outfile=trim(arg)
      i=i+2
    case ('-o2', '--out2')
      call get_command_argument(i+1, arg)
      cmdl%outfile_int=trim(arg)
      i=i+2
    case ('-c', '--cut')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%cutoff
      i=i+2
    case ('-ns', '--nspec')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%nspec
      i=i+2
    case ('-s', '--species')
      do j=1,cmdl%nspec
        call get_command_argument(i+j, arg)
        cmdl%speclist(j)=trim(arg)
      end do
      i=i+cmdl%nspec+1
    case ('-dr', '--binwidth')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%binwidth
      i=i+2
    case ('-sm', '--smooth')
      cmdl%smooth = .true.
      i=i+1
    case ('-H', '--noH')
      cmdl%noH = .true.
      i=i+1
    case ('-g', '--gwidth')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%gwidth
      i=i+2
    case ('-rm', '--rmin')
      call get_command_argument(i+1, arg)
      read(arg,*) cmdl%rmin
      i=i+2
    case ('--xsf')
      cmdl%xsf = .true.
      i=i+1
    case ('--outcar')
      cmdl%outcar = .true.
      i=i+1
    case ('--lammps')
      cmdl%lammps = .true.
      i=i+1
    case ('-f', '--freq')
      cmdl%freq = .true.
      i=i+1
    case default
      write(*,'(a,a,/)') 'Unrecognised command line option: ', arg
      call cmdl%print_help()
      stop
    end select
  end do

end subroutine get_args

end module cmdline
