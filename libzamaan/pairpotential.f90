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
  ! Stillinger-Weber variables
  real(double) :: A, B, p, q, lambda, gamma, swsigma, sweps, swcut, cos_th_ref

  contains
    procedure, public  :: read_pp_file
    procedure, public  :: init_pp
    procedure, public  :: init_lj
    procedure, public  :: init_morse
    procedure, public  :: init_sw
    procedure, public  :: pp_energy
    procedure, public  :: pp_force
    procedure, public  :: pp_force_and_energy
    procedure, public  :: lj_energy
    procedure, public  :: lj_force
    procedure, public  :: lj_force_and_energy
    procedure, public  :: morse_energy
    procedure, public  :: morse_force
    procedure, public  :: morse_force_and_energy
    procedure, private :: sw_v2
    procedure, private :: sw_f2
    procedure, private :: sw_h
    procedure, private :: sw_f3_e3
    procedure, public  :: sw_force_and_energy
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
    case ('stillinger-weber')
      read(funit,*) junk
      read(funit,*) pp%A
      read(funit,*) junk
      read(funit,*) pp%B
      read(funit,*) junk
      read(funit,*) pp%p
      read(funit,*) junk
      read(funit,*) pp%q
      read(funit,*) junk
      read(funit,*) pp%lambda
      read(funit,*) junk
      read(funit,*) pp%gamma
    end select

    close(funit)

  end subroutine read_pp_file

  subroutine init_pp(pp, filename, shift)

    ! passed variables
    class(type_pairpotential), intent(inout)  :: pp
    character(40), intent(in)                 :: filename
    logical, intent(in)                       :: shift

    ! local variables
    integer                                   :: i

    write(*,'(a)') "Initialising pair potential"

    call pp%read_pp_file(filename)

    write(*,'("Pair potential type    ",a)') pp%potential_type
    write(*,'("Pair potential cutoff  ",f8.4)') pp%r_cut
    write(*,'(2x,a)') "Sigma"
    do i=1,pp%ns
      write(*,'(4x,3f10.6)') pp%sigma(i,:)
    end do
    write(*,'(a)') "Epsilon"
    do i=1,pp%ns
      write(*,'(4x,3f10.6)') pp%epsilon(i,:)
    end do
    select case (pp%potential_type)
    case ("lennard-jones")
      call pp%init_lj(shift)
      write(*,'("Pair potential shift   ",l8)') shift
    case ("morse")
      call pp%init_morse(shift)
      write(*,'("Pair potential shift   ",l8)') shift
      write(*,'(2x,a)') "r_e"
      do i=1,pp%ns
        write(*,'(4x,3f10.6)') pp%r_e(i,:)
      end do
    case ('stillinger-weber')
      call pp%init_sw()
      write(*,'("A                      ",f8.4)') pp%A
      write(*,'("B                      ",f8.4)') pp%B
      write(*,'("p                      ",f8.4)') pp%p
      write(*,'("q                      ",f8.4)') pp%q
      write(*,'("lambda                 ",f8.4)') pp%lambda
      write(*,'("gamma                  ",f8.4)') pp%gamma
    end select
    write(*,*)

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

  ! Energy wrapper for arbitrary pair potential
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

  ! Force wrapper for arbitrary pair potential
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

  ! Force and energy wrapper for arbitrary pair potential
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

  ! Stillinger-Weber energy
  function sw_energy(sw, r_ij, s_i, s_j, shift) result(e)

    ! passed variables
    class(type_pairpotential), intent(in)  :: sw

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift
    real(double)                :: e

    ! local variables

  end function sw_energy

  ! Initialise Stillinger-Weber parameters
  subroutine init_sw(sw)

    ! passed variables
    class(type_pairpotential), intent(inout) :: sw

    sw%swsigma = sw%sigma(1,1)
    sw%sweps = sw%epsilon(1,1)
    sw%swcut = sw%r_cut(1,1)
    sw%cos_th_ref = -third

  end subroutine init_sw

  ! Stillinger-Weber 2-body potential energy
  function sw_v2(sw, rij, rijainv) result(v2)

    ! passed variables
    class(type_pairpotential), intent(in)  :: sw

    real(double), intent(in)    :: rij, rijainv
    real(double)                :: v2

    if (rij < sw%swcut) then
      v2 = sw%sweps*sw%A*(sw%B*rij**(-sw%p) - rij**(-sw%q))*exp(rijainv)
    else
      v2 = zero
    end if

  end function sw_v2

  ! Stillinger-Weber 2-body force
  function sw_f2(sw, rij, modrij, rijinv, rijainv, v2) result(f2)

    ! passed variables
    class(type_pairpotential), intent(in)   :: sw

    real(double), dimension(3), intent(in)  :: rij
    real(double), intent(in)                :: modrij, rijinv, rijainv, v2
    real(double), dimension(3)              :: f2

    if (modrij < sw%swcut) then
      f2 = v2*rij*rijinv*((sw%p*sw%B*modrij**(-sw%p-1) - &
                          sw%q*modrij**(-sw%q-1))/(sw%B*modrij**(-sw%p) - &
                          modrij**(-sw%q)) + rijainv**2)
    else
      f2 = zero
    end if

  end function sw_f2

  ! Stillinger-Weber three-body energy component
  function sw_h(sw, rinv1, rinv2, cos_th) result(h)

    ! passed variables
    class(type_pairpotential), intent(in)   :: sw
    real(double), intent(in)                :: rinv1, rinv2, cos_th
    real(double)                            :: h

    h = sw%sweps*sw%lambda*exp(sw%gamma*rinv1 + sw%gamma*rinv2) * &
        (cos_th-sw%cos_th_ref)**2

  end function sw_h

  ! Stillinger-Weber three-body force and energy
  subroutine sw_f3_e3(sw, rij, rjk, rik, rijinv, modrij, modrjk, modrik, &
                      rjkinv, rikinv, rijainv, rjkainv, rikainv, &
                      fi, fj, fk, h_jik, h_ijk, h_ikj)
    ! passed varialbes
    class(type_pairpotential), intent(in)   :: sw
    real(double), dimension(3), intent(in)  :: rij, rjk, rik
    real(double), intent(in)                :: modrij, modrjk, modrik, &
                                               rijinv, rjkinv, rikinv, &
                                               rijainv, rjkainv, rikainv
    real(double), dimension(3), intent(out) :: fi, fj, fk
    real(double), intent(out)               :: h_jik, h_ijk, h_ikj

    ! local variables
    real(double)                            :: cos_th_jik, cos_th_ijk, &
                                               cos_th_ikj

    ! jik permutation
    if ((modrij < sw%swcut) .and. (modrik < sw%swcut)) then
      cos_th_jik = dot_product(rij, rik)*rijinv*rikinv
      h_jik = sw%sw_h(rijainv, rikainv, cos_th_jik)
      fi = -sw%gamma*h_jik*(rij*rijinv*rijainv**2 + rik*rikinv*rikainv**2) + &
           two*sw%lambda*exp(sw%gamma*(rijainv + rikainv)) * &
           (cos_th_jik - sw%cos_th_ref) * &
           (rij*rijinv*rikinv + rik*rikinv*rijinv -  &
           (rij*rijinv*rikinv + rik*rikinv*rijinv)*cos_th_jik)
    else
      h_jik = zero
      fi = zero
    end if

    ! ijk permutation
    if ((modrij < sw%swcut) .and. (modrjk < sw%swcut)) then
      cos_th_ijk = dot_product(rij, rjk)*rijinv*rjkinv
      h_ijk = sw%sw_h(rijainv, rjkainv, cos_th_ijk)
      fj = -sw%gamma*h_ijk*(rij*rijinv*rijainv**2) + &
           two*sw%lambda*exp(sw%gamma*(rijainv + rjkainv)) * &
           (cos_th_ijk - sw%cos_th_ref) * &
           (-rjk*rjkinv*rijinv + rij*rijinv*rijinv*cos_th_ijk)
    else
      h_ijk = zero
      fj = zero
    end if

    ! ikj permutation
    if ((modrik < sw%swcut) .and. (modrjk < sw%swcut)) then
      cos_th_ikj = dot_product(rik, rjk)*rikinv*rjkinv
      h_ikj = sw%sw_h(rikainv, rjkainv, cos_th_ikj)
      fk = -sw%gamma*h_ikj*(rik*rikinv*rikainv**2) + &
           two*sw%lambda*exp(sw%gamma*(rikainv + rjkainv)) * &
           (cos_th_ikj - sw%cos_th_ref) * &
           (-rjk*rjkinv*rikinv + rik*rikinv*rikinv*cos_th_ikj)
    else
      h_ikj = zero
      fk = zero
    end if

  end subroutine sw_f3_e3


  ! Stillinger-Weber force and energy
  subroutine sw_force_and_energy(sw, rij, rik, rjk, ij_sw, jk_sw, ik_sw, &
                                 fi, fj, fk, v)

    ! passed variables
    class(type_pairpotential), intent(in)   :: sw
    real(double), dimension(3), intent(in)  :: rij, rjk, rik
    real(double), dimension(3), intent(out) :: fi, fj, fk
    real(double)                            :: v
    logical                                 :: ij_sw, jk_sw, ik_sw

    ! local variables
    real(double), dimension(3)              :: f3i, f3j, f3k, f2ij, f2jk, f2ik
    real(double)                            :: modrij, modrjk, modrik, &
                                               rijinv, rijainv, rjkinv, &
                                               rjkainv, rikinv, rikainv, &
                                               hjik, hijk, hikj, &
                                               v2ij, v2jk, v2ik, &
                                               cos_th_jik, cos_th_ijk, &
                                               cos_th_ikj

    v = zero
    fi = zero
    fj = zero
    fk = zero

    modrij = sqrt(sum(rij**2))
    modrjk = sqrt(sum(rjk**2))
    modrik = sqrt(sum(rik**2))
    if (modrij < sw%swcut) then
      rijinv = one/modrij
      rijainv = one/(modrij-sw%swcut)
      if (ij_sw) then
        v2ij = sw%sw_v2(modrij, rijainv)
        f2ij = sw%sw_f2(rij, modrij, rijinv, rijainv, v2ij)
      end if
    end if
    if (modrjk < sw%swcut) then
      rjkinv = one/modrjk
      rjkainv = one/(modrjk-sw%swcut)
      if (jk_sw) then
        v2jk = sw%sw_v2(modrjk, rjkainv)
        f2jk = sw%sw_f2(rjk, modrjk, rjkinv, rjkainv, v2jk)
      end if
    end if
    if (modrik < sw%swcut) then
      rikinv = one/modrik
      rikainv = one/(modrik-sw%swcut)
      if (ik_sw) then
        v2ik = sw%sw_v2(modrik, rikainv)
        f2ik = sw%sw_f2(rik, modrik, rikinv, rikainv, v2ik)
      end if
    end if

    call sw%sw_f3_e3(rij, rjk, rik, rijinv, modrij, modrjk, modrik, &
                     rjkinv, rikinv, rijainv, rjkainv, rikainv, &
                     f3i, f3j, f3k, hjik, hijk, hikj)

    fi = fi + f2ij + f2ik + f3i
    fj = fj - f2ij + f2jk + f3j
    fk = fk - f2ik - f2jk + f3k
    v = v + v2ij + v2jk + v2ik + hjik + hijk + hikj

  end subroutine sw_force_and_energy

end module pairpotential
