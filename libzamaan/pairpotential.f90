module pairpotential

use datatypes
use constants
use cell

implicit none

private
public :: type_pairpotential

type type_pairpotential
  character(40) :: potential_type
  integer       :: ns, npair, ntriple
  logical       :: threebody, shift
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
  real(double) :: A, B, p, q, lambda, gamma, swsigma, sigmainv, sweps, &
                  swcut, cos_th_ref

  contains
    procedure, public  :: read_pp_file
    procedure, public  :: init_pp
    procedure, public  :: init_lj
    procedure, public  :: init_morse
    procedure, public  :: init_sw
    procedure, public  :: pp_energy
    procedure, public  :: pp_force
    procedure, private :: pp_force_and_energy_component
    procedure, private :: pp_force_and_energy
    procedure, private :: mb_force_and_energy
    procedure, private :: lj_energy
    procedure, private :: lj_force
    procedure, private :: lj_force_and_energy
    procedure, private :: morse_energy
    procedure, private :: morse_force
    procedure, private :: morse_force_and_energy
    procedure, private :: sw_v2
    procedure, private :: sw_f2
    procedure, private :: sw_2body_force_and_energy
    procedure, private :: sw_h
    procedure, private :: sw_3body_force_and_energy
    procedure, private :: sw_force_and_energy
    procedure, public  :: get_force_and_energy
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

    pp%shift = shift
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
      call pp%init_lj
      write(*,'("Pair potential shift   ",l8)') pp%shift
    case ("morse")
      call pp%init_morse
      write(*,'("Pair potential shift   ",l8)') pp%shift
      write(*,'(2x,a)') "r_e"
      do i=1,pp%ns
        write(*,'(4x,3f10.6)') pp%r_e(i,:)
      end do
    case ('stillinger-weber')
      call pp%init_sw
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
  subroutine init_lj(pp)

    class(type_pairpotential), intent(inout)  :: pp

    integer                       :: i,j

    allocate(pp%s6(pp%ns,pp%ns), pp%s12(pp%ns,pp%ns), pp%c_f(pp%ns,pp%ns), &
                   pp%c_e(pp%ns,pp%ns))
    pp%s6 = pp%sigma**6
    pp%s12 = pp%s6**2
    pp%c_f = 24.0_double*pp%epsilon
    pp%c_e = 4.0_double*pp%epsilon

    ! If the potential is truncated there is a discontinuity at r_cut.
    ! Avoid by shifting the energy/force by adding U(r_cut) or F(r_cut)
    if (pp%shift .eqv. .true.) then
      allocate(pp%e_shift(pp%ns,pp%ns), pp%f_shift(pp%ns,pp%ns))
      pp%e_shift = zero
      pp%f_shift = zero
      do i=1,pp%ns
        do j=1,pp%ns
          pp%e_shift(i,j) = pp%lj_energy(pp%r_cut(i,j), i, j)
          pp%f_shift(i,j) = pp%lj_force(pp%r_cut(i,j), i, j)
        end do
      end do
    end if
  end subroutine init_lj

  ! Initialise a Morse potential
  subroutine init_morse(pp)

    class(type_pairpotential), intent(inout)  :: pp

    integer                       :: i,j

    allocate(pp%exp_re(pp%ns,pp%ns), pp%exp_re2(pp%ns,pp%ns), &
             pp%c_f(pp%ns,pp%ns))
    pp%exp_re = exp(pp%sigma*pp%r_e)
    pp%exp_re2 = pp%exp_re**2
    pp%c_f = two*pp%sigma*pp%epsilon

    ! If the potential is truncated there is a discontinuity at r_cut.
    ! Avoid by shifting the energy/force by adding U(r_cut) or F(r_cut)
    if (pp%shift .eqv. .true.) then
      allocate(pp%e_shift(pp%ns,pp%ns), pp%f_shift(pp%ns,pp%ns))
      pp%e_shift = zero
      pp%f_shift = zero
      do i=1,pp%ns
        do j=1,pp%ns
          pp%e_shift(i,j) = pp%morse_energy(pp%r_cut(i,j), i, j)
          pp%f_shift(i,j) = pp%morse_force(pp%r_cut(i,j), i, j)
        end do
      end do
    end if
  end subroutine init_morse

  ! Energy wrapper for arbitrary pair potential
  function pp_energy(pp, r_ij, s_i, s_j) result(e)

    ! passed variables
    class(type_pairpotential), intent(inout) :: pp

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: e

    select case (pp%potential_type)
    case ("lennard-jones")
      e = pp%lj_energy(r_ij, s_i, s_j)
    case ("morse")
      e = pp%morse_energy(r_ij, s_i, s_j)
    end select

  end function pp_energy

  ! Force wrapper for arbitrary pair potential
  function pp_force(pp, r_ij, s_i, s_j) result(f)

    ! passed variables
    class(type_pairpotential), intent(inout) :: pp

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: f

    select case (pp%potential_type)
    case ("lennard-jones")
      f = pp%lj_force(r_ij, s_i, s_j)
    case ("morse")
      f = pp%morse_force(r_ij, s_i, s_j)
    end select

  end function pp_force

  ! Force and energy contribution from one pair
  subroutine pp_force_and_energy_component(pp, r_ij, s_i, s_j, f, e)

    ! passed variables
    class(type_pairpotential), intent(inout) :: pp
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    real(double), intent(out)   :: f
    real(double), intent(out)   :: e

    select case (pp%potential_type)
    case ("lennard-jones")
      call pp%lj_force_and_energy(r_ij, s_i, s_j, f, e)
    case ("morse")
      call pp%morse_force_and_energy(r_ij, s_i, s_j, f, e)
    end select

  end subroutine pp_force_and_energy_component

  ! Total force and energy
  subroutine pp_force_and_energy(pp, p, v, f)

    ! passed variables
    class(type_pairpotential), intent(inout)  :: pp
    type(type_cell), intent(inout)            :: p
    real(double), intent(out)                 :: v
    real(double), dimension(:,:), intent(out) :: f

    ! local variables
    integer                                   :: i, j, s_i, s_j
    real(double), dimension(3)                :: rij, fij
    real(double)                              :: modrij, modf, vij

    call p%invert_lat
    call p%cell_cart2frac
    call p%get_vt

    v = zero
    f = zero

    do i=1,p%nat
      s_i = p%spec_int(i)
      do j=i+1,p%nat
        s_j = p%spec_int(j)
        fij = zero
        vij = zero

        if (p%dt(i,j) < pp%r_cut(s_i,s_j)) then
          call pp%pp_force_and_energy_component(p%dt(i,j), s_i, s_j, modf, vij)
          fij = modf*p%vt(i,j,:)/p%dt(i,j)
        end if
        f(i,:) = f(i,:) + fij
        f(j,:) = f(j,:) - fij
        v = v + vij

      end do
    end do

  end subroutine pp_force_and_energy

  ! Force and energy wrapper for manybody potential
  subroutine mb_force_and_energy(pp, p, v, f)

    ! passed variables
    class(type_pairpotential), intent(inout)  :: pp
    type(type_cell), intent(inout)            :: p
    real(double), intent(out)                 :: v
    real(double), dimension(:,:), intent(out) :: f

    select case(pp%potential_type)
    case('stillinger-weber')
      call pp%sw_force_and_energy(p, v, f)
    end select

  end subroutine mb_force_and_energy

  ! Lennard-Jones energy
  ! U = 4*epsil*((sigma/r)^12 - (sigma/r)^6)
  function lj_energy(lj, r_ij, s_i, s_j) result(e)

    ! passed variables
    class(type_pairpotential), intent(inout) :: lj

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: e
    ! local variables
    real(double)                :: inv_r6, inv_r12

    inv_r6 = one/(r_ij**6)
    inv_r12 = one/(r_ij**12)

    e = lj%c_e(s_i,s_j)*(lj%s12(s_i,s_j)*inv_r12-lj%s6(s_i,s_j)*inv_r6)
    if (lj%shift .eqv. .true.) then
        e = e - lj%e_shift(s_i,s_j)
    end if

  end function lj_energy

! Lennard-Jones force (derivative of energy)
! U = 24*epsil/sigma*(2*(sigma/r)^13 - (sigma/r)^7)
  function lj_force(lj, r_ij, s_i, s_j) result(f)

    ! passed variables
    class(type_pairpotential), intent(inout) :: lj

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: f
    ! local variables
    real(double)                :: inv_r7, inv_r13

    inv_r7 = one/(r_ij**7)
    inv_r13 = one/(r_ij**13)

    f = lj%c_f(s_i,s_j)*(2*lj%s12(s_i,s_j)*inv_r13-lj%s6(s_i,s_j)*inv_r7)
    if (lj%shift .eqv. .true.) then
      f = f - lj%f_shift(s_i,s_j)
    end if

  end function lj_force

  ! Compute force and energy in the same loop
  subroutine lj_force_and_energy(lj, r_ij, s_i, s_j, f, e)

    ! passed variables
    class(type_pairpotential), intent(inout) :: lj
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    real(double), intent(out)   :: f
    real(double), intent(out)   :: e

    ! local variables
    real(double)                :: inv_r, inv_r6

    inv_r = one/r_ij
    inv_r6 = inv_r**6

    e = lj%c_e(s_i,s_j)*lj%s6(s_i,s_j)*inv_r6*(inv_r6*lj%s6(s_i,s_j)-one)
    f = -lj%c_f(s_i,s_j)*inv_r*inv_r6*lj%s6(s_i,s_j)* &
        (two*lj%s6(s_i,s_j)*inv_r6 - one)
    if (lj%shift .eqv. .true.) then
      f = f - lj%f_shift(s_i,s_j)
      e = e - lj%e_shift(s_i,s_j)
    end if
  end subroutine lj_force_and_energy

  ! Morse energy
  ! U = eps[exp(-2*sigma(r-r_e))-2exp(-sigma(r-r_e))]
  function morse_energy(mo, r_ij, s_i, s_j) result(e)

    ! passed variables
    class(type_pairpotential), intent(inout) :: mo

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: e
    ! local variables
    real(double)                :: exp_1, exp_2

    exp_1 = exp(-mo%sigma(s_i,s_j)*r_ij)
    exp_2 =exp_1**2

    e = mo%epsilon(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - &
                             two*exp_1*mo%exp_re(s_i,s_j))

    if (mo%shift .eqv. .true.) then
        e = e - mo%e_shift(s_i,s_j)
    end if

  end function morse_energy

! Morse force (derivative of energy)
! U = 2*sigma*eps[exp(-2*sigma(r-r_e)) - exp(-sigma(r-r_e))]
  function morse_force(mo, r_ij, s_i, s_j) result(f)

    ! passed variables
    class(type_pairpotential), intent(inout) :: mo

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: f
    ! local variables
    real(double)                :: exp_1, exp_2

    exp_1 = exp(-mo%sigma(s_i,s_j)*r_ij)
    exp_2 =exp_1**2

    f = mo%c_f(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - exp_1*mo%exp_re(s_i,s_j))

    if (mo%shift .eqv. .true.) then
      f = f - mo%f_shift(s_i,s_j)
    end if

  end function morse_force

  ! Compute force and energy in the same loop
  subroutine morse_force_and_energy(mo, r_ij, s_i, s_j, f, e)

    ! passed variables
    class(type_pairpotential), intent(inout) :: mo
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    real(double), intent(out)   :: f
    real(double), intent(out)   :: e

    ! local variables
    real(double)                :: exp_1, exp_2

    exp_1 = exp(-mo%sigma(s_i,s_j)*r_ij)
    exp_2 =exp_1**2

    e = mo%epsilon(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - &
                             two*exp_1*mo%exp_re(s_i,s_j))
    f = mo%c_f(s_i,s_j)*(exp_2*mo%exp_re2(s_i,s_j) - exp_1*mo%exp_re(s_i,s_j))

    if (mo%shift .eqv. .true.) then
      f = f - mo%f_shift(s_i,s_j)
      e = e - mo%e_shift(s_i,s_j)
    end if
  end subroutine morse_force_and_energy

  ! Initialise Stillinger-Weber parameters
  subroutine init_sw(sw)

    ! passed variables
    class(type_pairpotential), intent(inout) :: sw

    sw%swsigma = sw%sigma(1,1)
    sw%sweps = sw%epsilon(1,1)
    sw%swcut = sw%r_cut(1,1)
    sw%cos_th_ref = -third
    sw%npair = 0
    sw%ntriple = 0
    sw%threebody = .true.

  end subroutine init_sw

  ! Stillinger-Weber 2-body potential energy
  function sw_v2(sw, rij, rijainv) result(v2)

    ! passed variables
    class(type_pairpotential), intent(inout) :: sw

    real(double), intent(in)    :: rij, rijainv
    real(double)                :: v2

    if (rij < sw%swcut) then
      v2 = sw%sweps*sw%A*(sw%B*rij**(-sw%p) - rij**(-sw%q))*exp(rijainv)
    end if

  end function sw_v2

  ! Stillinger-Weber 2-body force
  function sw_f2(sw, rij, modrij, rijinv, rijainv, v2) result(f2)

    ! passed variables
    class(type_pairpotential), intent(inout) :: sw

    real(double), dimension(3), intent(in)   :: rij
    real(double), intent(in)                 :: modrij, rijinv, rijainv, v2
    real(double), dimension(3)               :: f2

    sw%npair = sw%npair + 1
    if (modrij < sw%swcut) then
      f2 = -v2*rij*rijinv*(rijinv*(sw%p*sw%B*modrij**(-sw%p) - &
                          sw%q*modrij**(-sw%q))/(sw%B*modrij**(-sw%p) - &
                          modrij**(-sw%q)) + rijainv**2)
      ! f2 = sw%A*(sw%B*sw%p*modrij**(-sw%p) - &
        ! sw%q*modrij**(-sw%q))*rijinv*exp(rijainv) - &
        ! rijainv**2*sw%A*(sw%B*modrij**(-sw%p) - &
        ! modrij**(-sw%q))*exp(rijainv)
    end if

  end function sw_f2

  subroutine sw_2body_force_and_energy(sw, rij, fij, vij)

    ! passed variables
    class(type_pairpotential), intent(inout) :: sw
    real(double), intent(in), dimension(3)   :: rij
    real(double), dimension(3), intent(out)  :: fij
    real(double), intent(out)                :: vij

    ! local variables
    real(double)                             :: modrij, rijinv, rijainv

    vij = zero
    fij = zero

    modrij = sqrt(sum(rij**2))
    if (modrij < sw%swcut) then
      rijinv = one/modrij
      rijainv = one/(modrij-sw%swcut)
      vij = sw%sw_v2(modrij, rijainv)
      fij = sw%sw_f2(rij, modrij, rijinv, rijainv, vij)
    end if

  end subroutine sw_2body_force_and_energy

  ! Stillinger-Weber three-body energy component
  function sw_h(sw, rinv1, rinv2, cos_th) result(h)

    ! passed variables
    class(type_pairpotential), intent(inout)  :: sw
    real(double), intent(in)                :: rinv1, rinv2, cos_th
    real(double)                            :: h

    h = sw%sweps*sw%lambda*exp(sw%gamma*rinv1 + sw%gamma*rinv2) * &
        (cos_th-sw%cos_th_ref)**2

  end function sw_h

  ! Stillinger-Weber three-body force and energy
  subroutine sw_3body_force_and_energy(sw, rij, rjk, rik, modrij, modrjk, &
                      modrik, rijinv, rjkinv, rikinv, rijainv, rjkainv, &
                      rikainv, fi, fj, fk, h_jik, h_ijk, h_ikj)
    ! passed varialbes
    class(type_pairpotential), intent(inout) :: sw
    real(double), dimension(3), intent(in)   :: rij, rjk, rik
    real(double), intent(in)                 :: modrij, modrjk, modrik, &
                                                rijinv, rjkinv, rikinv, &
                                                rijainv, rjkainv, rikainv
    real(double), dimension(3), intent(out)  :: fi, fj, fk
    real(double), intent(out)                :: h_jik, h_ijk, h_ikj
 
    ! local variables
    real(double)                             :: cos_th_jik, cos_th_ijk, &
                                                cos_th_ikj

    sw%ntriple = sw%ntriple + 1

    ! jik permutation
    if ((modrij < sw%swcut) .and. (modrik < sw%swcut)) then
      cos_th_jik = dot_product(rij, rik)*rijinv*rikinv
      h_jik = sw%sw_h(rijainv, rikainv, cos_th_jik)
      fi = -sw%gamma*h_jik*(rij*rijinv*rijainv**2 + rik*rikinv*rikainv**2) + &
           two*sw%lambda*exp(sw%gamma*(rijainv + rikainv)) * &
           (cos_th_jik - sw%cos_th_ref) * &
           (rij*rijinv*rikinv + rik*rikinv*rijinv -  &
           (rij*rijinv*rikinv + rik*rikinv*rijinv)*cos_th_jik)
    end if

    ! ijk permutation
    if ((modrij < sw%swcut) .and. (modrjk < sw%swcut)) then
      cos_th_ijk = dot_product(rij, rjk)*rijinv*rjkinv
      h_ijk = sw%sw_h(rijainv, rjkainv, cos_th_ijk)
      fj = -sw%gamma*h_ijk*(rij*rijinv*rijainv**2) + &
           two*sw%lambda*exp(sw%gamma*(rijainv + rjkainv)) * &
           (cos_th_ijk - sw%cos_th_ref) * &
           (rjk*rjkinv*rijinv + rij*rijinv*rijinv*cos_th_ijk)
    end if

    ! ikj permutation
    if ((modrik < sw%swcut) .and. (modrjk < sw%swcut)) then
      cos_th_ikj = dot_product(rik, rjk)*rikinv*rjkinv
      h_ikj = sw%sw_h(rikainv, rjkainv, cos_th_ikj)
      fk = -sw%gamma*h_ikj*(rik*rikinv*rikainv**2) + &
           two*sw%lambda*exp(sw%gamma*(rikainv + rjkainv)) * &
           (cos_th_ikj - sw%cos_th_ref) * &
           (-rjk*rjkinv*rikinv + rik*rikinv*rikinv*cos_th_ikj)
    end if

  end subroutine sw_3body_force_and_energy

  subroutine sw_force_and_energy(sw, p, v, f)

    ! passed variables
    class(type_pairpotential), intent(inout)  :: sw
    type(type_cell), intent(inout)            :: p
    real(double), dimension(:,:), intent(out) :: f
    real(double), intent(out)                 :: v

    ! local variables
    integer                         :: i, j, k, s_i, s_j, s_k
    real(double), dimension(3)      :: fij, f3i, f3j, f3k
    real(double)                    :: vij,  hijk, hjik, hikj
    real(double)                    :: rijinv, rijainv, rikinv, rikainv, &
                                       rjkinv, rjkainv

    call p%invert_lat
    call p%cell_cart2frac
    call p%get_vt

    v = zero
    f = zero

    do i=1,p%nat
      s_i = p%spec_int(i)
      do j=i+1,p%nat
        s_j = p%spec_int(j)
        fij = zero
        vij = zero

        ! 2 body force and energy
        if (p%dt(i,j) < sw%swcut) then
          rijinv = one/p%dt(i,j)
          rijainv = one/(p%dt(i,j)-sw%swcut)
          vij = sw%sw_v2(p%dt(i,j), rijainv)
          fij = sw%sw_f2(p%vt(i,j,:), p%dt(i,j), rijinv, rijainv, vij)
        end if
        f(i,:) = f(i,:) + fij
        f(j,:) = f(j,:) - fij
        v = v + vij

        if (sw%threebody) then
          do k=j+1,p%nat
            s_k = p%spec_int(k)
            f3i = zero
            f3j = zero
            f3k = zero
            hijk = zero
            hikj = zero
            hjik = zero

            if (p%dt(j,k) < sw%swcut) then
              rjkinv = one/p%dt(j,k)
              rjkainv = one/(p%dt(j,k)-sw%swcut)
            end if

            if (p%dt(i,k) < sw%swcut) then
              rikinv = one/p%dt(i,k)
              rikainv = one/(p%dt(i,k)-sw%swcut)
            end if

            call sw%sw_3body_force_and_energy(p%vt(i,j,:), p%vt(j,k,:), &
              p%vt(i,k,:), p%dt(i,j), p%dt(j,k), p%dt(i,k), rijinv, &
              rjkinv, rikinv, rijainv, rjkainv, rikainv, f3i, f3j, f3k, &
              hjik, hijk, hikj)
            f(i,:) = f(i,:) + f3i
            f(j,:) = f(j,:) + f3j
            f(k,:) = f(k,:) + f3k
            v = v + hijk + hjik + hikj
          end do
        end if
      end do
    end do

  end subroutine sw_force_and_energy

  ! Wrapper for all potentials
  subroutine get_force_and_energy(pp, p, v, f)

    ! passed variables
    class(type_pairpotential), intent(inout)  :: pp
    type(type_cell), intent(inout)            :: p
    real(double), dimension(:,:), intent(out) :: f
    real(double), intent(out)                 :: v

    select case(pp%potential_type)
    case ("lennard-jones")
      call pp%pp_force_and_energy(p, v, f)
    case ("morse")
      call pp%pp_force_and_energy(p, v, f)
    case('stillinger-weber')
      call pp%mb_force_and_energy(p, v, f)
    case default
      stop "Unknown forcefield"
    end select

  end subroutine get_force_and_energy

end module pairpotential
