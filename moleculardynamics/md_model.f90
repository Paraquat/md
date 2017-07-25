module md_model

use datatypes

implicit none

type type_model
  
  ! Positions, velocities, forces
  real(double), dimension(:,:), allocatable :: r, r0
  real(double), dimension(:,:), allocatable :: rcart, rcart0
  real(double), dimension(:,:), allocatable :: v, v0
  real(double), dimension(:,:), allocatable :: f, f0

  ! Species
  integer, allocatable, dimension(:)        :: species
  character(2), allocatable, dimension(:)   :: spec_label
  real(double), allocatable, dimension(:)   :: mass

  ! box
  real(double), dimension(3,3)              :: h, h0

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
  real(double)                              :: T, T_ext
  real(double)                              :: P, P_ext
  real(double)                              :: vol, vol0
  real(double)                              :: Ek, Ek0
  real(double)                              :: Ep, Ep0
  real(double)                              :: enthalpy, enthalpy0
  real(double)                              :: pV, pV0
  real(double)                              :: H_prime
  real(double), dimension(3,3)              :: stress

  ! Thermostat
  character(40)                             :: thermo_type
  integer                                   :: n_nhc
  real(double), dimension(:), allocatable   :: Q
  real(double), dimension(:), allocatable   :: eta
  real(double), dimension(:), allocatable   :: v_eta
  real(double)                              :: G_eta
  real(double)                              :: E_nhc

  ! Barostat
  character(40)                             :: baro_type
  real(double)                              :: W_eps
  real(double)                              :: eps
  real(double)                              :: v_eps
  real(double)                              :: G_eps
  real(double)                              :: E_baro
  real(double), dimension(3,3)              :: c_g
  real(double)                              :: E_box

contains
  procedure :: init_model
  procedure :: dump_mdl_atom_arr
  procedure :: dump_frame
  procedure :: mdl_stat_dump

end type type_model

contains

  subroutine init_model(mdl, nat, nspec)

    ! passed variables
    class(type_model), intent(inout)  :: mdl
    integer, intent(in)               :: nat, nspec

    mdl%nat = nat
    mdl%nspec = nspec
    mdl%ndof = 3*nat

    allocate(mdl%r(nat,3), mdl%r0(nat,3))
    allocate(mdl%rcart(nat,3), mdl%rcart0(nat,3))
    allocate(mdl%v(nat,3), mdl%v0(nat,3))
    allocate(mdl%f(nat,3), mdl%f0(nat,3))
    allocate(mdl%species(nat))
    allocate(mdl%spec_label(nspec))
    allocate(mdl%mass(nspec))

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
  subroutine mdl_stat_dump(mdl, iunit, step)

    ! passed variables
    class(type_model), intent(inout)  :: mdl
    integer, intent(in)               :: iunit
    integer, intent(in)               :: step

    if (step == 0) then
      select case (mdl%ensemble)
      case ('nve')
        write(iunit,'(a10,5a16)') "step", "pe", "ke", "total", "T", "P"
      case ('nvt')
        if (mdl%thermo_type == 'nhc') then
          write(iunit,'(a10,6a16)') "step", "pe", "ke", "nhc", "total", "T", "P"
        else
          write(iunit,'(a10,5a16)') "step", "pe", "ke", "total", "T", "P"
        end if
      case ('npt')
        if (mdl%thermo_type == 'nhc') then
          write(iunit,'(a10,9a16)') "step", "pe", "ke", "nhc", "box", "pV", "total", "T", "P", "V"
        end if
      case ('nph')
        write(iunit,'(a10,8a16)') "step", "pe", "ke", "box", "pV", "total", "T", "P", "V"
    end select
    end if
    select case (mdl%ensemble)
    case ('nve')
      write(iunit,'(i10,5e16.6)') step, mdl%Ep, mdl%Ek, mdl%H_prime, &
                                  mdl%T, mdl%P
    case ('nvt')
      if (mdl%thermo_type == 'nhc') then
        write(iunit,'(i10,6e16.6)') step, mdl%Ep, mdl%Ek, mdl%E_nhc, &
                                    mdl%H_prime, mdl%T, mdl%P
      else
        write(iunit,'(i10,5e16.6)') step, mdl%Ep, mdl%Ek, mdl%H_prime, &
                                    mdl%T, mdl%P
      end if
    case ('npt')
      if (mdl%thermo_type == 'nhc') then
        write(iunit,'(i10,9e16.6)') step, mdl%Ep, mdl%Ek, mdl%E_nhc, &
                                    mdl%E_box, mdl%pV, mdl%H_prime, &
                                    mdl%T, mdl%P, mdl%V
      end if
    case ('nph')
        write(iunit,'(i10,8e16.6)') step, mdl%Ep, mdl%Ek, mdl%E_box, &
                                    mdl%pV, mdl%H_prime, mdl%T, mdl%P, mdl%V
    end select

  end subroutine mdl_stat_dump


end module md_model