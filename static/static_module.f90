module static_module

use datatypes
use constants
use vector
use cell
use pairpotential
use rng

implicit none

type type_static
  type(type_cell) :: p
  type(type_pp)   :: pp
  integer         :: nspec, nat, ndof, iprint
  logical         :: shift
  character(40)   :: position_file, dump_file
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:,:)   :: mass
  real(double), allocatable, dimension(:,:)   :: f
  real(double), dimension(3,3)                :: P_int, P_ext
  real(double)   :: pe

  contains
    procedure :: init_static
    procedure :: get_force_and_energy
    procedure :: get_pressure
    procedure :: dump_atom_arr
end type type_static

contains

  ! Initialise variables/velocities/allocate matrices for MD run
  subroutine init_static(static, init_cell, pp, init_cell_cart, shift)

    ! passed variables
    class(type_static), intent(inout)   :: static
    type(type_cell), intent(in)     :: init_cell
    type(type_pp), intent(in)       :: pp
    logical, intent(in)             :: init_cell_cart
    logical, intent(in)             :: shift

    ! local variables
    integer                         :: i, j

    static%position_file = "trajectory.xsf"
    static%dump_file = "dump.out"
    write(*,'(a)') "Starting static run"
    write(*,*)
    static%p = init_cell
    if (init_cell_cart .eqv. .false.) call static%p%cell_frac2cart
    static%nspec = static%p%nspec
    static%pp = pp
    static%iprint = 0
    static%nat = init_cell%nat
    static%shift = shift

    allocate(static%species(static%nat))
    allocate(static%f(static%nat,3))
    static%species = init_cell%spec_int

    write(*,'(a)') "Species and pair potential details:"
    write(*,'("Pair potential cutoff  ",f8.4)') static%pp%r_cut
    write(*,'("Pair potential shift   ",l8)') static%shift
    write(*,'(2x,a)') "Sigma"
    do i=1,static%nspec
      write(*,'(4x,3f10.6)') static%pp%sigma(i,:)
    end do
    write(*,'(a)') "Epsilon"
    do i=1,static%nspec
      write(*,'(4x,3f10.6)') static%pp%epsilon(i,:)
    end do
    write(*,*)
    write(*,'(2x,a)') "label, species, count"
    do i=1,static%nspec
      write(*,'(4x,i4,a4,i8)') static%p%spec_int(i), static%p%spec(i), &
                               static%p%spec_count(i)
    end do
    write(*,*)
    write(*,'(2x,a)') "Initial unit cell:"
    do i=1,3
      write(*,'(4x,3f12.6)') static%p%h(i,:)
    end do
    write(*,*)
    write(*,'(2x,a)') "Initial atomic positions:"
    do i=1,static%nat
      write(*,'(4x,2i6,3f14.8)') i, static%species(i), static%p%rcart(i,:)
    end do

  end subroutine init_static

  ! Compute the potential energy and force on each atom
  subroutine get_force_and_energy(static)

    ! passed variables
    class(type_static), intent(inout)   :: static

    ! local variables
    integer                         :: iat, jat, s_i, s_j
    real(double), dimension(3)      :: r_ij_cart
    real(double)                    :: mod_r_ij, mod_f, pe

    write(*,*)
    write(*,'(a)') "Computing forces and energy"
    static%f = zero
    static%pe = zero

    do iat=1,static%nat
      do jat=1,static%nat
        if (iat == jat) cycle
        s_i = static%species(iat)
        s_j = static%species(jat)
        r_ij_cart = static%p%mic(static%p%rcart(iat,:), &
                                 static%p%rcart(jat,:))
        mod_r_ij = modulus(r_ij_cart)
        if (mod_r_ij < static%pp%r_cut(s_i,s_j)) then
          r_ij_cart = norm(r_ij_cart)
          call static%pp%lj_force_and_energy(mod_r_ij, s_i, s_j, static%shift, &
                                          mod_f, pe)
          static%f(iat,:) = static%f(iat,:) + mod_f*r_ij_cart
          static%pe = static%pe + pe
        end if
      end do
    end do
    static%pe = static%pe*half
    write(*,*)
    write(*,'("  Potential energy = ",e16.8)') static%pe
    write(*,*)
    do iat=1,static%nat
      write(*,'("  Force: ",i8,3e20.10)') iat, static%f(iat,:)
    end do
  end subroutine get_force_and_energy

  ! Compute the pressure using the virial
  subroutine get_pressure(static)

    ! passed variables
    class(type_static), intent(inout)   :: static

    ! local variables
    integer                         :: i

    do i=1,static%nat
    end do

  end subroutine get_pressure

  ! Dump the an atom array (position, force, velocity)
  subroutine dump_atom_arr(static, iou, arr)

    ! passed variables
    class(type_static), intent(inout)             :: static
    integer, intent(in)                       :: iou
    real(double), dimension(:,:), intent(in)  :: arr

    ! local variables
    integer                         :: i
    character(80)                   :: fmt

    fmt = "(2i5,3e20.10)"
    do i=1,static%nat
      write(iou,fmt) i, static%species(i), arr(i,:)
    end do

  end subroutine dump_atom_arr

end module static_module
