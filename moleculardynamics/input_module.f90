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
    logical         :: restart = .false.    ! restart from checkpoint
    character(3)    :: ensemble = 'nve'
    real(double)    :: dt                   ! time step
    integer         :: natoms
    integer         :: nsteps
    integer         :: dump_freq            ! time steps between dumps
    integer         :: iprint = 0           ! print level
    integer         :: cp_freq = 100        ! time steps between checkpoint
    logical         :: comv = .true.        ! remove centre of mass velocity
    logical         :: cart = .false.       ! input structure in Cartesian coord
    character(40)   :: v_distr = 'maxwell-boltzmann' ! velocity distribution
    character(40)   :: pbc_method = 'frac'
    logical         :: rdf = .false.        ! compute g(r)?
    real(double)    :: rdfcut = 10.0        ! g(r) cutoff
    real(double)    :: dr                   ! g(r) bin width
    real(double)    :: rmin                 ! g(r) minimum
    real(double)    :: gwidth               ! g(r) Gaussian smoothing width

    ! Pair potential paramters
    logical         :: shift = .false.

    ! Thermostat parameters
    character(40)   :: thermo_type = 'none'
    real(double)    :: T_ext                ! external temperature
    real(double)    :: tau_T = one          ! thermostat time period
    integer         :: n_nhc = 5            ! length of Nose-Hoover chain
    integer         :: n_mts = 1            ! multiple time step order
    integer         :: n_ys = 1             ! Yoshida-Suzuki order
    logical         :: cell_nhc = .true.    ! use separate NHC for cell?
    real(double), dimension(:), allocatable :: nhc_mass
    real(double), dimension(:), allocatable :: cell_nhc_mass

    ! Barostat paramters
    character(40)   :: baro_type = 'none'
    real(double)    :: P_ext                ! external pressure
    real(double)    :: box_mass = one       ! box mass for extended-Lagrangian
    real(double)    :: tau_P = one          ! barostat time period
    real(double)    :: bulkmod = 100.0_double ! bulk modulus estimate

  contains
    procedure, public :: read_input

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
      write(*,'("thermo_type         ",a20)') inp%thermo_type
    case ('T_ext')
      read(param,*) inp%T_ext
      write(*,'("T_ext               ",f20.8)') inp%T_ext
    case ('tau_T')
      read(param,*) inp%tau_T
      write(*,'("tau_T               ",f20.8)') inp%tau_T
    case ('n_nhc')
      read(param,*) inp%n_nhc
      write(*,'("n_nhc               ",i20)') inp%n_nhc
      allocate(inp%nhc_mass(inp%n_nhc))
      allocate(inp%cell_nhc_mass(inp%n_nhc))
      inp%nhc_mass = one
      inp%cell_nhc_mass = one
    case ('v_distr')
      read(param,*) inp%v_distr
      write(*,'("v_distr             ",a20)') inp%v_distr
    case('cell_nhc')
      read(param,*) inp%cell_nhc
      write(*,'("v_distr             ",l20)') inp%cell_nhc
    case ('nhc_mass')
      if (inp%n_nhc > 0) then
        read(param,*) inp%nhc_mass
        write(*,'("nhc_mass          ",f20.6)') inp%nhc_mass
      end if
    case ('cell_nhc_mass')
      if (inp%n_nhc > 0) then
        read(param,*) inp%cell_nhc_mass
        write(*,'("cell_nhc_mass      ",f20.6)') inp%cell_nhc_mass
      end if
    case ('baro_type')
      read(param,*) inp%baro_type
      write(*,'("baro_type           ",a20)') inp%baro_type
    case ('P_ext')
      read(param,*) inp%P_ext
      write(*,'("P_ext               ",f20.8)') inp%P_ext
    case ('box_mass')
      read(param,*) inp%box_mass
      write(*,'("box_mass            ",f20.8)') inp%box_mass
    case ('tau_P')
      read(param,*) inp%tau_P
      write(*,'("tau_P               ",f20.8)') inp%tau_P
    case ('restart')
      read(param,*) inp%restart
      write(*,'("restart             ",l20)') inp%restart
    case('checkpoint')
      read(param,*) inp%cp_freq
      write(*,'("checkpoint          ",i20)') inp%cp_freq
    case ('pbc_method')
      read(param,*) inp%pbc_method
      write(*,'("pbc_method          ",a20)') inp%pbc_method
    case ('bulkmod_est')
      read(param,*) inp%bulkmod
      write(*,'("bulkmod_est         ",f20.6)') inp%bulkmod
    case ('iprint')
      read(param,*) inp%iprint
      write(*,'("iprint              ",i20)') inp%iprint
    case ('rdf')
      read(param,*) inp%rdf
      write(*,'("rdf                 ",l20)') inp%rdf
    case ('rdfcut')
      read(param,*) inp%rdfcut
      write(*,'("rdfcut              ",f20.6)') inp%rdfcut
    case ('dr')
      read(param,*) inp%dr
      write(*,'("dr                  ",f20.6)') inp%dr
    case ('rmin')
      read(param,*) inp%rmin
      write(*,'("rmin                ",f20.6)') inp%rmin
    case ('gwidth')
      read(param,*) inp%gwidth
      write(*,'("gwidth              ",f20.6)') inp%gwidth
    end select

  end do
  write(*,*)

end subroutine read_input

end module input_module
