module input_module

  use datatypes
  use constants

  implicit none

  type md_input
    ! Files
    character(40)   :: sfile = 'cell.in'    ! Structure file
    character(40)   :: mdfile = 'md.in'     ! MD parameters file
    character(40)   :: ppfile = 'pp.in'     ! pair potential file
    character(40)   :: dumpfile = 'dump.out'
    character(40)   :: statfile = 'stat.out'

    ! General simulation parameters
    character(3)    :: ensemble = 'nve'
    real(double)    :: dt                   ! time step
    integer         :: natoms
    integer         :: nsteps
    integer         :: dump_freq            ! time steps between dumps
    logical         :: comv = .true.        ! remove centre of mass velocity
    logical         :: cart = .false.       ! input structure in Cartesian coord
    character(40)   :: v_distr = 'uniform'  ! Initial velocity distribution

    ! Pair potential paramters
    logical         :: shift = .false.

    ! Thermostat parameters
    character(40)   :: thermo_type = 'velocity_rescale'
    real(double)    :: T_ext                ! external temperature
    integer         :: tau_T = 1            ! thermostat time period
    integer         :: n_nhc = 5            ! lenght of Nose-Hoover chain

    ! Barostat paramters
    real(double)        :: P_ext

  contains
    procedure :: read_input

  end type md_input

contains

subroutine read_input(inp, filename)

  ! passed variables
  class(md_input), intent(inout)  :: inp
  character(40), intent(in)       :: filename

  ! local variables
  character(80)                   :: buffer, label, param
  integer                         :: ios, paramcount, fh, line, pos

  fh = 15
  open(fh, file=filename, status='old', iostat=ios)
  if (ios > 0) stop "Could not find md parameters file (md.in)"
  write(*,'("Reading MD parameters from ",a)') filename
  do
    read(fh, '(a)', iostat=ios) buffer
    if (ios > 0) then
      write(*,'(a)') 'Error reading input file'
      stop
    else if (ios < 0) then ! end of file
      exit
    end if

    ! Split the buffer on the first instance of whitespace
    pos=scan(buffer, ' ')
    label = buffer(:pos)
    param = trim(buffer(pos+1:))

    select case (label)
    case ('sfile')
      read(param,*) inp%sfile
      write(*,'("sfile               ",20a)') inp%sfile
    case ('dumpfile')
      read(param,*) inp%dumpfile
      write(*,'("dumpfile            ",20a)') inp%dumpfile
    case ('ppfile')
      read(param,*) inp%dumpfile
      write(*,'("ppfile              ",20a)') inp%ppfile
    case ('statfile')
      read(param,*) inp%statfile
      write(*,'("statfile            ",20a)') inp%statfile
    case ('ensemble')
      read(param,*) inp%ensemble
      write(*,'("ensemble            ",20a)') inp%ensemble
    case ('dt')
      read(param,*) inp%dt
      write(*,'("dt                  ",f20.8)') inp%dt
    case ('nsteps')
      read(param,*) inp%nsteps
      write(*,'("nsteps              ",i20)') inp%nsteps
    case ('natoms')
      read(param,*) inp%natoms
      write(*,'("natoms              ",i20)') inp%natoms
    case ('dump_freq')
      read(param,*) inp%dump_freq
      write(*,'("dump_freq           ",i20)') inp%dump_freq
    case ('comv')
      read(param,*) inp%comv
      write(*,'("dump_freq           ",l20)') inp%comv
    case ('cart')
      read(param,*) inp%cart
      write(*,'("cart                ",l20)') inp%cart
    case ('shift')
      read(param,*) inp%shift
      write(*,'("shift               ",l20)') inp%shift
    case ('thermo_type')
      read(param,*) inp%thermo_type
      write(*,'("thermo_type         ",l20)') inp%thermo_type
    case ('T_ext')
      read(param,*) inp%T_ext
      write(*,'("T_ext               ",l20)') inp%T_ext
    case ('tau_T')
      read(param,*) inp%tau_T
      write(*,'("tau_T               ",l20)') inp%tau_T
    case ('n_nhc')
      read(param,*) inp%n_nhc
      write(*,'("n_nhc               ",i20)') inp%n_nhc
    case ('P_ext')
      read(param,*) inp%P_ext
      write(*,'("P_ext               ",l20)') inp%P_ext
    case ('v_distr')
      read(param,*) inp%v_distr
      write(*,'("v_distr             ",a20)') inp%v_distr
    end select

  end do
  write(*,*)
end subroutine read_input

end module input_module
