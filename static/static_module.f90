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
  integer         :: nspec, nat, ndof, iprint, eos_max, ncalc, dos_ngrids
  logical         :: shift, rdf
  character(40)   :: position_file, dump_file, gr_file, runtype
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:)     :: modf, volume, energy
  real(double), allocatable, dimension(:,:)   :: mass
  real(double), allocatable, dimension(:,:)   :: f
  real(double), allocatable, dimension(:,:)   :: fc_matrix, dyn_matrix
  real(double), dimension(3,3)                :: static_stress, h_ref
  real(double)   :: pe, eos_incr, displ, r_fc

  contains
    procedure :: init_static
    procedure :: get_force_and_energy
    procedure :: run_eos
    procedure :: run_fdphonon
    procedure :: get_fc_matrix
    procedure :: truncate_fc_matrix
    procedure :: acoustic_sum_rule
    procedure :: get_phonons
    procedure :: get_phonon_dos
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
    case ('fdphonon')
      allocate(sp%fc_matrix(sp%nat*3,sp%nat*3))
      sp%fc_matrix = zero
      sp%displ = inp%displ
      sp%r_fc = inp%r_fc
      sp%dos_ngrids = inp%dos_ngrids
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

  ! Compute the equaion of state
  subroutine run_eos(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    character(40)                       :: eosfile
    integer                             :: icalc, eosunit
    real(double), dimension(:), allocatable :: sf

    eosfile = 'eos.dat'
    eosunit = 21
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
      call sp%p%invert_lat
      call sp%p%cell_frac2cart
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

  ! Compute phonons
  subroutine run_fdphonon(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    real(double), dimension(:,:), allocatable :: phonon_modes
    real(double), dimension(:), allocatable   :: omega_sq, phonon_dos_total
    integer                                   :: pdos_unit
    character(40)                             :: pdos_file
    real(double)                              :: bin_mid, grid_sp
    integer                                   :: i

    if (sp%iprint > 1) write(*,'(2x,a)') "run_fdphonon"

    allocate(omega_sq(3*sp%nat))
    allocate(phonon_modes(3*sp%nat,3*sp%nat))
    allocate(phonon_dos_total(sp%dos_ngrids))
    pdos_unit = 22
    pdos_file = "pdos.dat"

    call sp%get_fc_matrix
    call sp%truncate_fc_matrix
    call sp%acoustic_sum_rule
    call sp%get_phonons(phonon_modes, omega_sq)
    call sp%get_phonon_dos(omega_sq, phonon_dos_total, grid_sp)

    open(unit=pdos_unit,file=pdos_file,status='replace')
    do i=1,sp%dos_ngrids
      bin_mid = real(((i-1)/2), double)
      write(pdos_unit,'(2f20.12)') bin_mid*grid_sp, phonon_dos_total(i)
    end do
    close(pdos_unit)

  end subroutine run_fdphonon

  subroutine get_fc_matrix(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    integer                             :: i, j, k, l, i_fc, j_fc
    real(double), dimension(3)          :: pos_temp

    if (sp%iprint > 1) write(*,'(2x,a)') "get_fc_matrix"

    do i=1,sp%nat
      pos_temp = sp%p%rcart(i,:) 
      do j=1,3
        i_fc = 3*(i-1)+j
        sp%p%rcart(i,j) = sp%p%rcart(i,j) + sp%displ  
        write(*,'(2x,"Force calculation ",i6," of ",i6)') i_fc, 9*sp%nat**2
        call sp%get_force_and_energy
        do k=1,sp%nat
          do l=1,3
            j_fc = 3*(k-1)+l
            sp%fc_matrix(i_fc,j_fc) = sp%f(k,l)
          end do
        end do

        sp%p%rcart(i,:) = pos_temp
      end do
    end do
    write(*,'(2x,a)') "Phonon force calculations complete"
  end subroutine get_fc_matrix

  ! Truncate force constant matrix using fc range r_fc
  subroutine truncate_fc_matrix(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    integer                             :: iat, jat, i, j, idir, jdir

    if (sp%iprint > 1) write(*,'(2x,a)') "truncate_fc_matrix"

    do iat=1,sp%nat
      do jat=1,sp%nat
        if (sp%p%dt(iat,jat) > sp%r_fc) then
          i=3*iat-2
          j=3*jat-2
          do idir=0,2
            do jdir=0,2
              sp%fc_matrix(i+idir,j+jdir) = zero
              sp%fc_matrix(j+jdir,i+idir) = zero
            end do
          end do
        end if
      end do
    end do
  end subroutine truncate_fc_matrix

  ! Impose the acoustic sum rule
  subroutine acoustic_sum_rule(sp)

    ! passed variables
    class(type_static), intent(inout)   :: sp

    ! local variables
    integer                             :: iat, jat, idir, jdir, idm, jdm
    real(double)                        :: total
    real(double), dimension(3)          :: f_asr

    if (sp%iprint > 1) write(*,'(2x,a)') "acoustic_sum_rule"

    write(*,'(2x,a)') "Checking acoustic sum rule:"
    write(*,'(4x,a8,a3,a3,a12)') "atom", "i", "j", "fsum"
    do iat=1,sp%nat
      do idir=1,3
        do jdir=1,3
          idm = 3*(iat-1)+idir
          total = sum(sp%fc_matrix(idm,jdir::3))
          write(*,'(4x,i8,i3,i3,e16.8)') iat, idir, jdir, total
        end do
      end do
    end do

    ! Correction
    do idm=1,3*sp%nat
      do idir=1,3
        f_asr(idir) = sum(sp%fc_matrix(idm,idir::3))
      end do
      do jat=1,sp%nat
        do jdir=1,3
          jdm = 3*(jat-1)+jdir
          sp%fc_matrix(idm,jdm) = sp%fc_matrix(idm,jdm)-f_asr(jdir)/real(sp%nat,double)
        end do
      end do
    end do
    
    write(*,*)
    write(*,'(2x,a)') "After correction:"
    write(*,'(4x,a8,a3,a3,a12)') "atom", "i", "j", "fsum"
    do iat=1,sp%nat
      do idir=1,3
        do jdir=1,3
          idm = 3*(iat-1)+idir
          total = sum(sp%fc_matrix(idm,jdir::3))
          write(*,'(4x,i8,i3,i3,e16.8)') iat, idir, jdir, total
        end do
      end do
    end do

  end subroutine acoustic_sum_rule

  ! Get the dynamical matrix and diagonalise it
  subroutine get_phonons(sp, dyn_mat_sym, omega_sq)

    ! passed variables
    class(type_static), intent(inout)         :: sp
    real(double), dimension(:,:), intent(out) :: dyn_mat_sym
    real(double), dimension(:), intent(out)   :: omega_sq

    ! local variables
    integer                                   :: dir_i, dir_j, lwork, info, &
                                                 spec_i, spec_j, i, j
    real(double)                              :: E0, dm_const, mass_i, mass_j
    real(double), dimension(:), allocatable   :: work

    if (sp%iprint > 1) write(*,'(2x,a)') "get_phonons"

    lwork=9*sp%nat
    allocate(work(lwork))
    work = zero

    ! Construct dynamical matrix in real space
    do dir_i=1,3*sp%nat
      i = ((dir_i-1)/3) + 1 
      spec_i = sp%p%spec_int(i)
      mass_i = sp%p%mass(spec_i)
      do dir_j=1,3*sp%nat
        j = ((dir_j-1)/3) + 1 
        spec_j = sp%p%spec_int(j)
        mass_j = sp%p%mass(spec_j)
        dm_const = one/sqrt(mass_i*mass_j)/sp%displ
        sp%fc_matrix(dir_i,dir_j) = sp%fc_matrix(dir_i,dir_j)*dm_const
        dyn_mat_sym(dir_j,dir_i) = sp%fc_matrix(dir_i,dir_j)
      end do
    end do
  
    ! symmetrise the dynamical matrix
    E0 = zero
    do dir_i=1,3*sp%nat
      do dir_j=1,3*sp%nat
        E0 = E0 + (dyn_mat_sym(dir_j,dir_i)-sp%fc_matrix(dir_j,dir_i))**2
      end do
    end do
    write(*,'(2x,"Norm of anti-symm is ",f20.12)') sqrt(E0)
    dyn_mat_sym = half*(dyn_mat_sym+sp%fc_matrix)

    E0 = zero
    do dir_i=1,3*sp%nat
      do dir_j=1,3*sp%nat
        E0 = E0 + (dyn_mat_sym(dir_j,dir_i)+sp%fc_matrix(dir_j,dir_i))**2
      end do
    end do
    write(*,'(2x,"Norm of symm is      ",f20.12)') sqrt(E0)

    ! Diagonalise the dynamical matrix
    write(*,'(2x,a)') "Diagonalising the dynamical matrix"
    call dsyev('V', 'U', 3*sp%nat, dyn_mat_sym, 3*sp%nat, omega_sq, work, &
               lwork, info)
    write(*,'(2x,a)') "Finished diagonalising"
    write(*,'(4x,"Info: ",i4)') info
    write(*,'(4x,a)') "Eigenvalues:"
    do dir_i=1,3*sp%nat
      write(*,'(6x,i4,f13.7)') dir_i, omega_sq(dir_i)
    end do

  end subroutine get_phonons

  ! Get the phonon DOS
  subroutine get_phonon_dos(sp, omega_sq, phonon_dos_total, dos_grid_sp)

    ! passed variables
    class(type_static), intent(inout)       :: sp
    real(double), dimension(:), intent(in)  :: omega_sq
    real(double), dimension(:), intent(out) :: phonon_dos_total
    real(double), intent(out)               :: dos_grid_sp

    ! local variables
    integer                                 :: neig, nim, i, idos
    real(double)                            :: emax, eig, norm_fac

    if (sp%iprint > 1) write(*,'(2x,a)') "get_phonon_dos"

    neig = 3*sp%nat
    nim = 0

    emax = sqrt(maxval(omega_sq))
    dos_grid_sp = emax/real(sp%dos_ngrids, double)
    phonon_dos_total = zero

    write(*,'(2x,"DOS grid spacing:  ",f12.4)') dos_grid_sp
    write(*,'(2x,"Number of grids:   ",i12)') sp%dos_ngrids
    write(*,'(2x,"Maximum frequency: ",f12.4)') emax

    do i=1,neig
      if (omega_sq(i) > zero) then
        eig = sqrt(omega_sq(i))
        idos = int(eig/dos_grid_sp)
        if (idos > 0) then
          phonon_dos_total(idos) = phonon_dos_total(idos) + 1
        end if
      else if (omega_sq(i) < -small) then
        nim = nim + 1
      end if
    end do

    ! Normalisation: integral = 1/ndof = neig
    norm_fac = one/real(neig,double)
    phonon_dos_total = phonon_dos_total*norm_fac

    if (nim > 1) then
      write(*,'(2x,"Warning: ",i4," imaginary frequencies")') nim
    end if

  end subroutine get_phonon_dos

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
    case('fdphonon')
      write(*,'(a)') "Calculating finite difference phonons"
      if (sp%rdf) then
        call sp%pd%update_rdist(sp%p)
        call sp%pd%norm_rdist
        call sp%pd%write_gr(sp%gr_file)
      end if
      call sp%run_fdphonon
    end select
  end subroutine run_static

end module static_module