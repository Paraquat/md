module pairpotential

use datatypes
use constants

implicit none

private
public :: type_pairpotential

type type_pairpotential
  character(40) :: potential_type
  integer       :: ns
  real(double), allocatable, dimension(:,:)  :: r_cut
  character(2), allocatable, dimension(:)    :: species_label
  real(double), allocatable, dimension(:,:)  :: sigma
  real(double), allocatable, dimension(:,:)  :: epsilon
  real(double), allocatable, dimension(:,:)  :: r_e
  real(double), allocatable, dimension(:,:)  :: f_shift, e_shift
  ! Lennard-Jones variables
  real(double), allocatable, dimension(:,:)  :: s6,s12,c_f,c_e
  ! Morse variables
  real(double), allocatable, dimension(:,:)  :: exp_re, exp_re2

  contains
    procedure :: read_pp_file
    procedure :: init_pp
    procedure :: init_lj
    procedure :: init_morse
    procedure :: pp_energy
    procedure :: pp_force
    procedure :: pp_force_and_energy
    procedure :: lj_energy
    procedure :: lj_force
    procedure :: lj_force_and_energy
    procedure :: morse_energy
    procedure :: morse_force
    procedure :: morse_force_and_energy
end type type_pairpotential

contains


  ! Read a potential file and initialise
  ! Format:
  ! 1 number of atoms
  ! 2 list of species
  ! 3 Comment (# Cutoff)
  ! 4 r_cut matrix
  ! 5 Comment (# Epsilon)
  ! 6 epsilon matrix
  ! 7 Comment (# Sigma)
  ! 8 sigma matrix
  ! 9 Comment (# r_e)
  ! 10 r_e matrix
  subroutine read_pp_file(pp, filename)

    class(type_pairpotential), intent(inout)  :: pp

    character(40), intent(in)     :: filename

    integer                       :: funit, ios, i
    character(80)                 :: junk

    funit = 103

    open(funit, file=filename, iostat=ios)
    read(funit,*) pp%potential_type
    pp%potential_type = trim(pp%potential_type)
    read(funit,*) pp%ns

    allocate(pp%sigma(pp%ns,pp%ns), pp%epsilon(pp%ns,pp%ns), &
             pp%r_cut(pp%ns,pp%ns), pp%species_label(pp%ns))
    select case(pp%potential_type)
    case ('morse')
      allocate(pp%r_e(pp%ns,pp%ns))
    end select

    read(funit,*) pp%species_label

    read(funit,*) junk
    do i=1,pp%ns
      read(funit,*) pp%r_cut(i,:)
    end do

    read(funit,*) junk
    do i=1,pp%ns
      read(funit,*) pp%epsilon(i,:)
    end do

    read(funit,*) junk
    do i=1,pp%ns
      read(funit,*) pp%sigma(i,:)
    end do

    select case(pp%potential_type)
    case ('morse')
      read(funit,*) junk
      do i=1,pp%ns
        read(funit,*) pp%r_e(i,:)
      end do
    end select

    close(funit)

  end subroutine read_pp_file

  subroutine init_pp(pp, filename, shift)

    class(type_pairpotential), intent(inout)  :: pp
    character(40), intent(in)                 :: filename
    logical, intent(in)                       :: shift

    call pp%read_pp_file(filename)

    select case (pp%potential_type)
    case ("lennard-jones")
      call pp%init_lj(shift)
    case ("morse")
      call pp%init_morse(shift)
    end select

  end subroutine init_pp

  ! Initialise a Lennard-Jones potential
  subroutine init_lj(pp, shift)

    class(type_pairpotential), intent(inout)  :: pp
    logical, intent(in)           :: shift

    integer                       :: i,j

    allocate(pp%s6(pp%ns,pp%ns), pp%s12(pp%ns,pp%ns), pp%c_f(pp%ns,pp%ns), &
                   pp%c_e(pp%ns,pp%ns))
    pp%s6 = pp%sigma**6
    pp%s12 = pp%s6**2
    pp%c_f = 24.0_double*pp%epsilon
    pp%c_e = 4.0_double*pp%epsilon

    ! If the potential is truncated there is a discontinuity at r_cut.
    ! Avoid by shifting the energy/force by adding U(r_cut) or F(r_cut)
    if (shift .eqv. .true.) then
      allocate(pp%e_shift(pp%ns,pp%ns), pp%f_shift(pp%ns,pp%ns))
      pp%e_shift = zero
      pp%f_shift = zero
      do i=1,pp%ns
        do j=1,pp%ns
          pp%e_shift(i,j) = pp%lj_energy(pp%r_cut(i,j), i, j, .false.)
          pp%f_shift(i,j) = pp%lj_force(pp%r_cut(i,j), i, j, .false.)
        end do
      end do
    end if
  end subroutine init_lj

  ! Initialise a Morse potential
  subroutine init_morse(mo, shift)

    class(type_pairpotential), intent(inout)  :: mo
    logical, intent(in)           :: shift

    integer                       :: i,j

    allocate(mo%exp_re(mo%ns,mo%ns), mo%exp_re2(mo%ns,mo%ns), &
             mo%c_f(mo%ns,mo%ns))
    mo%exp_re = exp(mo%sigma*mo%r_e)
    mo%exp_re2 = mo%exp_re**2
    mo%c_f = two*mo%sigma*mo%epsilon

    ! If the potential is truncated there is a discontinuity at r_cut.
    ! Avoid by shifting the energy/force by adding U(r_cut) or F(r_cut)
    if (shift .eqv. .true.) then
      allocate(mo%e_shift(mo%ns,mo%ns), mo%f_shift(mo%ns,mo%ns))
      mo%e_shift = zero
      mo%f_shift = zero
      do i=1,mo%ns
        do j=1,mo%ns
          mo%e_shift(i,j) = mo%morse_energy(mo%r_cut(i,j), i, j, .false.)
          mo%f_shift(i,j) = mo%morse_force(mo%r_cut(i,j), i, j, .false.)
        end do
      end do
    end if
  end subroutine init_morse

  ! Lennard-Jones energy
  ! U = 4*epsil*((sigma/r)^12 - (sigma/r)^6)
  function pp_energy(pp, r_ij, s_i, s_j, shift) result(e)

    ! passed variables
    class(type_pairpotential), intent(in)  :: pp

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: e

    select case (pp%potential_type)
    case ("lennard-jones")
      e = pp%lj_energy(r_ij, s_i, s_j, shift)
    case ("morse")
      e = pp%morse_energy(r_ij, s_i, s_j, shift)
    end select

  end function pp_energy

! Lennard-Jones force (derivative of energy)
! U = 24*epsil/sigma*(2*(sigma/r)^13 - (sigma/r)^7)
  function pp_force(pp, r_ij, s_i, s_j, shift) result(f)

    ! passed variables
    class(type_pairpotential), intent(in)  :: pp

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: f

    select case (pp%potential_type)
    case ("lennard-jones")
      f = pp%lj_force(r_ij, s_i, s_j, shift)
    case ("morse")
      f = pp%morse_force(r_ij, s_i, s_j, shift)
    end select

  end function pp_force

  ! Compute force and energy in the same loop
  subroutine pp_force_and_energy(pp, r_ij, s_i, s_j, shift, f, e)

    ! passed variables
    class(type_pairpotential), intent(in)  :: pp
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift
    real(double), intent(out)   :: f
    real(double), intent(out)   :: e

    select case (pp%potential_type)
    case ("lennard-jones")
      call pp%lj_force_and_energy(r_ij, s_i, s_j, shift, f, e)
    case ("morse")
      call pp%morse_force_and_energy(r_ij, s_i, s_j, shift, f, e)
    end select

  end subroutine pp_force_and_energy

  ! Lennard-Jones energy
  ! U = 4*epsil*((sigma/r)^12 - (sigma/r)^6)
  function lj_energy(lj, r_ij, s_i, s_j, shift) result(e)

    ! passed variables
    class(type_pairpotential), intent(in)  :: lj

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: e
    ! local variables
    real(double)                :: inv_r6, inv_r12

    inv_r6 = one/(r_ij**6)
    inv_r12 = one/(r_ij**12)

    e = lj%c_e(s_i,s_j)*(lj%s12(s_i,s_j)*inv_r12-lj%s6(s_i,s_j)*inv_r6)
    if (shift .eqv. .true.) then
        e = e - lj%e_shift(s_i,s_j)
    end if

  end function lj_energy

! Lennard-Jones force (derivative of energy)
! U = 24*epsil/sigma*(2*(sigma/r)^13 - (sigma/r)^7)
  function lj_force(lj, r_ij, s_i, s_j, shift) result(f)

    ! passed variables
    class(type_pairpotential), intent(in)  :: lj

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: f
    ! local variables
    real(double)                :: inv_r7, inv_r13

    inv_r7 = one/(r_ij**7)
    inv_r13 = one/(r_ij**13)

    f = lj%c_f(s_i,s_j)*(2*lj%s12(s_i,s_j)*inv_r13-lj%s6(s_i,s_j)*inv_r7)
    if (shift .eqv. .true.) then
      f = f - lj%f_shift(s_i,s_j)
    end if

  end function lj_force

  ! Compute force and energy in the same loop
  subroutine lj_force_and_energy(lj, r_ij, s_i, s_j, shift, f, e)

    ! passed variables
    class(type_pairpotential), intent(in)  :: lj
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift
    real(double), intent(out)   :: f
    real(double), intent(out)   :: e

    ! local variables
    real(double)                :: inv_r, inv_r6

    inv_r = one/r_ij
    inv_r6 = inv_r**6

    e = lj%c_e(s_i,s_j)*lj%s6(s_i,s_j)*inv_r6*(inv_r6*lj%s6(s_i,s_j)-one)
    f = -lj%c_f(s_i,s_j)*inv_r*inv_r6*lj%s6(s_i,s_j)* &
        (two*lj%s6(s_i,s_j)*inv_r6 - one)
    if (shift .eqv. .true.) then
      f = f - lj%f_shift(s_i,s_j)
      e = e - lj%e_shift(s_i,s_j)
    end if
  end subroutine lj_force_and_energy

  ! Morse energy
  ! U = eps[exp(-2*sigma(r-r_e))-2exp(-sigma(r-r_e))]
  function morse_energy(mo, r_ij, s_i, s_j, shift) result(e)

    ! passed variables
    class(type_pairpotential), intent(in)  :: mo

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: e
    ! local variables
    real(double)                :: exp_1, exp_2

    exp_1 = exp(-mo%sigma(s_i,s_j)*r_ij)
    exp_2 =exp_1**2

    e = mo%epsilon(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - &
                             two*exp_1*mo%exp_re(s_i,s_j))

    if (shift .eqv. .true.) then
        e = e - mo%e_shift(s_i,s_j)
    end if

  end function morse_energy

! Morse force (derivative of energy)
! U = 2*sigma*eps[exp(-2*sigma(r-r_e)) - exp(-sigma(r-r_e))]
  function morse_force(mo, r_ij, s_i, s_j, shift) result(f)

    ! passed variables
    class(type_pairpotential), intent(in)  :: mo

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: f
    ! local variables
    real(double)                :: exp_1, exp_2

    exp_1 = exp(-mo%sigma(s_i,s_j)*r_ij)
    exp_2 =exp_1**2

    f = mo%c_f(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - exp_1*mo%exp_re(s_i,s_j))

    if (shift .eqv. .true.) then
      f = f - mo%f_shift(s_i,s_j)
    end if

  end function morse_force

  ! Compute force and energy in the same loop
  subroutine morse_force_and_energy(mo, r_ij, s_i, s_j, shift, f, e)

    ! passed variables
    class(type_pairpotential), intent(in)  :: mo
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift
    real(double), intent(out)   :: f
    real(double), intent(out)   :: e

    ! local variables
    real(double)                :: exp_1, exp_2

    exp_1 = exp(-mo%sigma(s_i,s_j)*r_ij)
    exp_2 =exp_1**2

    e = mo%epsilon(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - &
                             two*exp_1*mo%exp_re(s_i,s_j))
    f = mo%c_f(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - exp_1*mo%exp_re(s_i,s_j))

    if (shift .eqv. .true.) then
      f = f - mo%f_shift(s_i,s_j)
      e = e - mo%e_shift(s_i,s_j)
    end if
  end subroutine morse_force_and_energy

end module pairpotential
