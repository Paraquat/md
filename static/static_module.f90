module static_module

use datatypes
use constants
use vector
use cell
use pairpotential
use rng
use pairdist

implicit none

type type_static
  type(type_cell) :: p
  type(type_pairpotential)   :: pp
  type(type_pairdist) :: pd
  integer         :: nspec, nat, ndof, iprint, eos_max, ncalc
  logical         :: shift, rdf
  character(40)   :: position_file, dump_file, gr_file, runtype
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:)     :: modf, volume, energy
  real(double), allocatable, dimension(:,:)   :: mass
  real(double), allocatable, dimension(:,:)   :: f
  real(double), dimension(3,3)                :: static_stress, h_ref
  real(double)   :: pe, eos_incr

  contains
    procedure :: init_static
    procedure :: get_force_and_energy
    procedure :: run_eos
    procedure :: get_static_stress
    procedure :: dump_atom_arr
    procedure :: run_static
end type type_static

contains

  ! Initialise variables/velocities/allocate matrices for MD run
  subroutine init_static(sp, inp, init_cell, pp)

    use input_module, only: static_input

    ! passed variables
    class(type_static), intent(inout)     :: sp
    type(static_input), intent(in)        :: inp
    type(type_cell), intent(in)           :: init_cell
    type(type_pairpotential), intent(in)  :: pp

    ! local variables
    integer                         :: i, j

    sp%runtype = inp%runtype
    sp%position_file = "trajectory.xsf"
    sp%dump_file = "dump.out"
    sp%gr_file = 'gr.dat'
    write(*,'(a)') "Starting static run"
    write(*,*)
    sp%p = init_cell
    if (inp%cart .eqv. .false.) then
      call sp%p%cell_frac2cart
    else
      call sp%p%cell_cart2frac
    end if
    sp%nspec = sp%p%nspec
    sp%pp = pp
    sp%iprint = 0
    sp%nat = init_cell%nat
    sp%shift = inp%shift

    sp%rdf = inp%rdf
    if (inp%rdf) then
      call sp%pd%init_pd(sp%p, inp%rdfcut, inp%gwidth, inp%dr, inp%rmin)
      sp%pd%smooth_on = .false.
    end if

    allocate(sp%species(sp%nat))
    allocate(sp%f(sp%nat,3))
    allocate(sp%modf(sp%nat))
    sp%species = init_cell%spec_int

    write(*,*)
    write(*,'(2x,a)') "label, species, count"
    do i=1,sp%nspec
      write(*,'(4x,i4,a4,i8)') sp%p%spec_int(i), sp%p%spec(i), &
                               sp%p%spec_count(i)
    end do
    write(*,*)
    write(*,'(2x,a)') "Initial unit cell:"
    do i=1,3
      write(*,'(4x,3f12.6)') sp%p%h(i,:)
    end do
    write(*,*)
    write(*,'(2x,a)') "Initial atomic positions (fractional):"
    do i=1,sp%nat
      write(*,'(4x,2i6,3f14.8)') i, sp%species(i), sp%p%r(i,:)
    end do

    select case(sp%runtype)
    case ('static')
    case ('eos')
      sp%eos_incr = inp%eos_incr
      sp%eos_max = inp%eos_max
      sp%ncalc = 2*sp%eos_max + 1
      allocate(sp%volume(sp%ncalc))
      allocate(sp%energy(sp%ncalc))
    end select

  end subroutine init_static

  ! Compute the potential energy and force on each atom
  subroutine get_force_and_energy(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    integer                             :: iat, imax

    if (sp%iprint > 1) write(*,'(2x,a)') "get_force_and_energy"
    call sp%pp%get_force_and_energy(sp%p, sp%pe, sp%f)
    call sp%get_static_stress
    write(*,*)
    do iat=1,sp%nat
      write(*,'("  Force: ",i8,3e20.10)') iat, sp%f(iat,:)
      sp%modf(iat) = sqrt(sum(sp%f(iat,:)**2))
    end do
    imax = maxloc(sp%modf,1)
    write(*,*)
    write(*,'("  Maximum force on atom ", i6, ": ", 3f14.8)') &
      imax, sp%f(imax,:)
    write(*,'("  Potential energy = ",e16.8)') sp%pe
  end subroutine get_force_and_energy

  subroutine run_eos(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    character(40)                       :: eosfile
    integer                             :: icalc, eosunit
    real(double), dimension(:), allocatable :: sf

    eosfile = 'eos.dat'
    if (sp%iprint > 1) write(*,'(2x,a)') "run_eos"

    allocate(sf(sp%ncalc))
    sf(1) = one - real(sp%eos_max, double)*sp%eos_incr
    do icalc = 2,sp%ncalc
      sf(icalc) = sf(icalc-1) + sp%eos_incr
    end do
    sp%h_ref = sp%p%h
    call sp%p%invert_lat
    call sp%p%cell_frac2cart

    do icalc = 1,sp%ncalc
      write(*,'(2x,"  EOS energy calculation ",i4, " of ",i4)') &
        icalc, sp%ncalc
      sp%p%h = sp%h_ref*sf(icalc)
      call sp%get_force_and_energy
      sp%volume(icalc) = sp%p%volume()
      sp%energy(icalc) = sp%pe
      write(*,'(2x,"Volume: ",f12.4," Energy: ",e12.6)') &
        sp%volume(icalc), sp%energy(icalc)
    end do

    deallocate(sf)
    open(unit=eosunit,file=eosfile,status='replace')
    do icalc=1,sp%ncalc
      write(eosunit,'(2f20.10)') sp%volume(icalc), sp%energy(icalc)
    end do
    close(eosunit)
    write(*,'(2x,a)') "Equation of state calculation complete"

  end subroutine run_eos

  ! Compute the pressure using the virial
  subroutine get_static_stress(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    integer                         :: i

    if (sp%iprint > 1) write(*,'(2x,a)') "get_static_stress"
    sp%static_stress = zero

    do i=1,sp%nat
      sp%static_stress = sp%static_stress + &
                          tensor_product(sp%f(i,:), sp%p%rcart(i,:))
    end do
    sp%static_stress = sp%static_stress/sp%p%volume()

    write(*,'(2x,a)') "Static stress:"
    do i=1,3
      write(*,'(4x,3f14.8)') sp%static_stress(i,:)
    end do

  end subroutine get_static_stress

  ! Dump the an atom array (position, force, velocity)
  subroutine dump_atom_arr(sp, iou, arr)

    ! passed variables
    class(type_static), intent(inout)         :: sp
    integer, intent(in)                       :: iou
    real(double), dimension(:,:), intent(in)  :: arr

    ! local variables
    integer                         :: i
    character(80)                   :: fmt

    fmt = "(2i5,3e20.10)"
    do i=1,sp%nat
      write(iou,fmt) i, sp%species(i), arr(i,:)
    end do

  end subroutine dump_atom_arr

  ! Wrapper for all static calculation types
  subroutine run_static(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    if (sp%iprint > 1) write(*,'(2x,a)') "run_static"

    select case (sp%runtype)
    case('static')
      write(*,'(a)') "Performing static energy calculation"
      call sp%get_force_and_energy
      if (sp%rdf) then
        call sp%pd%update_rdist(sp%p)
        call sp%pd%norm_rdist
        call sp%pd%write_gr(sp%gr_file)
      end if
    case('eos')
      write(*,'(a)') "Calculating equation of state"
      call sp%run_eos
    end select
  end subroutine run_static

end module static_module