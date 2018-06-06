module md_model

use datatypes
use input_module
use cell
use md_control

implicit none

type type_model

  logical                                   :: append
  
  ! Positions, velocities, forces
  real(double), dimension(:,:), pointer     :: r
  real(double), dimension(:,:), pointer     :: rcart
  real(double), dimension(:,:), pointer     :: v
  real(double), dimension(:,:), pointer     :: f

  ! real(double), dimension(:,:), allocatable :: r, r0
  ! real(double), dimension(:,:), allocatable :: rcart, rcart0
  ! real(double), dimension(:,:), allocatable :: v, v0
  ! real(double), dimension(:,:), allocatable :: f, f0

  ! Species
  integer, dimension(:), allocatable        :: species
  character(2), dimension(:), pointer       :: spec_label
  real(double), dimension(:), pointer       :: mass

  ! box
  real(double), dimension(:,:), pointer     :: h, h0

  !  MD variables
  integer                                   :: step
  integer                                   :: nstep
  real(double)                              :: dt
  character(3)                              :: ensemble

  ! Microscopic variables
  integer                                   :: nat
  integer                                   :: nspec
  integer                                   :: ndof

  ! Thermodynamic variables
  real(double), pointer                     :: T_int, T_ext
  real(double), pointer                     :: P_int, P_ext
  real(double), pointer                     :: vol, vol_ref
  real(double)                              :: E_k
  real(double)                              :: E_p
  real(double)                              :: PV
  real(double)                              :: H_prime
  real(double), dimension(:,:), pointer     :: kinetic_stress, static_stress
  real(double), dimension(:,:), pointer     :: stress

  ! Thermostat
  type(type_thermostat), pointer            :: th
  character(40), pointer                    :: thermo_type
  integer, pointer                          :: n_nhc
  real(double), dimension(:), pointer       :: eta
  real(double), dimension(:), pointer       :: eta_cell
  real(double), dimension(:), pointer       :: v_eta
  real(double), dimension(:), pointer       :: v_eta_cell
  real(double), dimension(:), pointer       :: m_nhc
  real(double), dimension(:), pointer       :: m_nhc_cell
  real(double), dimension(:), pointer       :: G_eta
  real(double), dimension(:), pointer       :: G_eta_cell
  real(double), pointer                     :: E_nhc

  ! Barostat
  character(40), pointer                    :: baro_type
  real(double), pointer                     :: m_box
  real(double), pointer                     :: eps
  real(double), pointer                     :: v_eps
  real(double), pointer                     :: G_eps
  real(double), pointer                     :: E_baro
  real(double), dimension(:), pointer       :: Q
  real(double), dimension(:), pointer       :: v_Q
  real(double), dimension(:), pointer       :: G_Q
  real(double), dimension(:,:), pointer     :: v_h
  real(double), dimension(:,:), pointer     :: G_h
  real(double), pointer                     :: E_box

contains
  procedure :: init_model
  procedure :: dump_mdl_atom_arr
  procedure :: dump_frame
  procedure :: stat_dump

end type type_model

contains

  subroutine init_model(mdl, inp, cell, v, f, th, baro, ndof)

    ! passed variables
    class(type_model), intent(inout)                  :: mdl
    type(md_input), intent(in)                        :: inp
    type(type_cell), intent(in), target               :: cell
    real(double), dimension(:,:), intent(in), target  :: v, f
    type(type_thermostat), intent(in), target         :: th
    type(type_barostat), intent(in), target           :: baro
    integer, intent(in)                               :: ndof

    if (inp%restart) then
      mdl%append = .true.
    else
      mdl%append = .false.
    end if

    mdl%ndof = ndof
    mdl%nat = inp%natoms
    mdl%dt = inp%dt
    mdl%ensemble = inp%ensemble

    mdl%E_nhc => th%e_nhc
    mdl%E_box => baro%ke_box

    allocate(mdl%species(mdl%nat))
    mdl%species = cell%spec_int
    mdl%r => cell%r
    mdl%rcart => cell%rcart
    mdl%h => cell%h
    mdl%v => v
    mdl%f => f

    mdl%T_int       => th%T_int
    mdl%T_ext       => th%T_ext
    mdl%thermo_type => th%th_type
    mdl%n_nhc       => th%n_nhc
    mdl%eta         => th%eta
    mdl%v_eta       => th%v_eta
    mdl%G_eta       => th%G_nhc
    mdl%m_nhc       => th%m_nhc
    if (th%cell_nhc) then
      mdl%eta_cell    => th%eta_cell
      mdl%v_eta_cell  => th%v_eta_cell
      mdl%G_eta_cell  => th%G_nhc_cell
      mdl%m_nhc_cell  => th%m_nhc_cell
    end if

    mdl%P_int       => baro%P_int
    mdl%P_ext       => baro%P_ext
    mdl%vol         => baro%V
    mdl%vol_ref     => baro%V_ref
    mdl%stress      => baro%stress
    mdl%kinetic_stress => baro%kinetic_stress
    mdl%static_stress  => baro%static_stress
    mdl%baro_type   => baro%baro_type
    mdl%m_box       => baro%m_box
    mdl%eps         => baro%eps
    mdl%v_eps       => baro%v_eps
    mdl%G_eps       => baro%G_eps
    mdl%Q           => baro%Q
    mdl%v_Q         => baro%v_Q
    mdl%G_Q         => baro%G_Q
    mdl%v_h         => baro%v_h
    mdl%G_h         => baro%G_h

  end subroutine init_model

 ! Dump the an atom array (position, force, velocity)
  subroutine dump_mdl_atom_arr(mdl, iou, arr)

    ! passed variables
    class(type_model), intent(inout)          :: mdl
    integer, intent(in)                       :: iou
    real(double), dimension(:,:), intent(in)  :: arr

    ! local variables
    integer                           :: i
    character(80)                     :: fmt

    fmt = "(2i5,3e20.10)"
    do i=1,mdl%nat
      write(iou,fmt) i, mdl%species(i), arr(i,:)
    end do

  end subroutine dump_mdl_atom_arr

  subroutine dump_frame(mdl, iunit, step)

    ! passed variables
    class(type_model), intent(inout)  :: mdl
    integer, intent(in)               :: iunit
    integer, intent(in)               :: step

    ! local variables
    integer                           :: i

    write(iunit,'("frame ",i8)') step
    write(iunit,'(a)') "cell_vectors"
    do i=1,3
      write(iunit,'(3f12.6)') mdl%h(i,:)
    end do
    write(iunit,'(a)') "end cell_vectors"
    write(iunit,'(a)') "stress_tensor"
    do i=1,3
      write(iunit,'(3f12.6)') mdl%stress(i,:)
    end do
    write(iunit,'(a)') "end stress_tensor"
    write(iunit,'(a)') "positions"
    call mdl%dump_mdl_atom_arr(iunit, mdl%rcart)
    write(iunit,'(a)') "end positions"
    write(iunit,'(a)') "velocities"
    call mdl%dump_mdl_atom_arr(iunit, mdl%v)
    write(iunit,'(a)') "end velocities"
    write(iunit,'(a)') "forces"
    call mdl%dump_mdl_atom_arr(iunit, mdl%f)
    write(iunit,'(a)') "end forces"
    write(iunit,'(a)') "end frame"

  end subroutine dump_frame

 ! Dump thermodyanmic statistics
  subroutine stat_dump(mdl, iunit, step)

    ! passed variables
    class(type_model), intent(inout)  :: mdl
    integer, intent(in)               :: iunit
    integer, intent(in)               :: step

    if (step == 0) then
      select case (mdl%ensemble)
      case ('nve')
        write(iunit,'(a10,5a16)') "step", "pe", "ke", "H'", "T", "P"
      case ('nvt')
        if (mdl%thermo_type == 'nhc') then
          write(iunit,'(a10,6a16)') "step", "pe", "ke", "nhc", "H'", "T", "P"
        else
          write(iunit,'(a10,5a16)') "step", "pe", "ke", "H'", "T", "P"
        end if
      case ('npt')
        if (mdl%thermo_type == 'nhc') then
          write(iunit,'(a10,9a16)') "step", "pe", "ke", "nhc", "box", "pV", "H'", "T", "P", "V"
        end if
      case ('nph')
        write(iunit,'(a10,8a16)') "step", "pe", "ke", "box", "pV", "H'", "T", "P", "V"
    end select
    end if
    select case (mdl%ensemble)
    case ('nve')
      write(iunit,'(i10,5e16.6)') step, mdl%E_p, mdl%E_k, mdl%H_prime, &
                                  mdl%T_int, mdl%P_int
    case ('nvt')
      if (mdl%thermo_type == 'nhc') then
        write(iunit,'(i10,6e16.6)') step, mdl%E_p, mdl%E_k, mdl%E_nhc, &
                                    mdl%H_prime, mdl%T_int, mdl%P_int
      else
        write(iunit,'(i10,5e16.6)') step, mdl%E_p, mdl%E_k, mdl%H_prime, &
                                    mdl%T_int, mdl%P_int
      end if
    case ('npt')
      if (mdl%thermo_type == 'nhc') then
        write(iunit,'(i10,9e16.6)') step, mdl%E_p, mdl%E_k, mdl%E_nhc, &
                                    mdl%E_box, mdl%pV, mdl%H_prime, &
                                    mdl%T_int, mdl%P_int, mdl%vol
      end if
    case ('nph')
        write(iunit,'(i10,8e16.6)') step, mdl%E_p, mdl%E_k, mdl%E_box, &
                                    mdl%pV, mdl%H_prime, mdl%T_int, &
                                    mdl%P_int, mdl%vol
    end select

  end subroutine stat_dump


end module md_model