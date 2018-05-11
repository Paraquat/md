module md_control

use datatypes
use constants
use rng
use maths,        only: poly_sinhx_x, mmult3x3, mmult3x3_loop
use input_module, only: md_input

implicit none

type type_thermostat

  real(double)        :: T_int, T_ext, k_B_md, ke_atoms
  integer             :: nat, ndof, iprint
  integer             :: n_ys    ! Yoshida-Suzuki order
  integer             :: n_nhc   ! number of nh chains
  integer             :: n_mts   ! number of NHC loops per time step
  real(double)        :: e_nhc   ! NHC contribution to the conserved term
  real(double)        :: tau_T   ! temperature coupling time period
  real(double)        :: dt      ! time step
  real(double)        :: lambda  ! Berendsen scaling factor
  character(40)       :: th_type
  logical             :: cell_nhc ! whether to use a separate NHC for the cell
  character(3)        :: ensemble
  real(double), dimension(:), allocatable   :: eta   ! thermostat k position
  real(double), dimension(:), allocatable   :: eta_cell
  real(double), dimension(:), allocatable   :: v_eta ! thermostat k velocity
  real(double), dimension(:), allocatable   :: v_eta_cell
  real(double), dimension(:), allocatable   :: G_nhc ! force on thermostat k
  real(double), dimension(:), allocatable   :: G_nhc_cell
  real(double), dimension(:), allocatable   :: m_nhc ! thermostat k mass
  real(double), dimension(:), allocatable   :: m_nhc_cell
  real(double), dimension(:), allocatable   :: dt_ys ! Y-S time steps

contains

  procedure :: init_thermostat
  procedure :: init_ys
  procedure :: get_temperature
  procedure :: velocity_rescale
  procedure :: get_berendsen_thermo_sf
  procedure :: berendsen_thermo_propagate
  procedure :: init_nhc
  procedure :: update_G_nhc
  procedure :: propagate_eta
  procedure :: propagate_v_eta_lin
  procedure :: propagate_v_eta_exp
  procedure :: propagate_nvt_nhc
  procedure :: get_nhc_ke
  procedure :: dump_thermo_state

end type type_thermostat

type type_barostat

  character(40)                 :: baro_type
  character(3)                  :: ensemble
  real(double)                  :: dt
  integer                       :: iprint
  integer                       :: nat
  integer                       :: ndof
  real(double)                  :: k_B_md
  real(double)                  :: bulkmod

  real(double)                  :: P_int    ! internal pressure
  real(double)                  :: P_ext    ! applied pressure
  real(double)                  :: V        ! volume
  real(double)                  :: V_ref    ! reference volume
  real(double)                  :: odnf
  real(double), dimension(3,3)  :: h        ! current cell
  real(double), dimension(3,3)  :: h_0      ! refernce cell
  real(double), dimension(3,3)  :: h_scale  ! cells scaling factors
  real(double), dimension(3,3)  :: stress   ! stress tensor
  real(double), dimension(3,3)  :: stress_ext
  real(double), dimension(3,3)  :: kinetic_stress ! ke contribution to stress tensor
  real(double), dimension(3,3)  :: static_stress ! static contribution to stress tensor

  ! MTTK variables

  real(double)                  :: m_box      ! box mass
  real(double)                  :: eps
  real(double)                  :: v_eps
  real(double)                  :: G_eps
  real(double), dimension(3)    :: Q
  real(double), dimension(3)    :: v_Q
  real(double), dimension(3)    :: G_Q
  real(double), dimension(3,3)  :: v_h      ! box velocity
  real(double), dimension(3,3)  :: G_h      ! box force
  real(double), dimension(3)    :: G_nhc_1  !
  real(double), dimension(3,3)  :: ident    ! 3x3 identity matrix
  real(double), dimension(3,3)  :: onfm     ! ident*(1/ndof)
  real(double), dimension(3,3)  :: c_g      ! eigenvectors 
  real(double), dimension(3,3)  :: I_e
  real(double), dimension(3,3)  :: I_s
  real(double), dimension(3)    :: lambda   ! eigenvalues
  real(double)                  :: ke_box
  real(double)                  :: ke_atoms

  ! Berendsen variables
  real(double)                  :: tau_P    ! pressure coupling time constatn
  real(double)                  :: mu       ! Berendsen scaling factor
  real(double), dimension(3,3)  :: mu_h     ! Berendsen scaling matrix

contains

  procedure, public  :: init_barostat
  procedure, public  :: get_berendsen_baro_sf
  procedure, public  :: berendsen_baro_propagate
  procedure, public  :: get_box_ke
  procedure, public  :: update_stress
  procedure, public  :: get_stress_and_pressure
  procedure, public  :: update_G_box
  procedure, public  :: propagate_eps_lin
  procedure, private :: propagate_v_eps_lin
  procedure, private :: propagate_v_eps_exp
  procedure, public  :: propagate_box
  procedure, private :: diag_vbox
  procedure, public  :: get_Ie
  procedure, public  :: get_Is
  procedure, private :: propagate_v_h_1
  procedure, private :: propagate_v_h_2
  procedure, public  :: propagate_v_ions
  procedure, public  :: propagate_r_ions
  procedure, public  :: propagate_h
  procedure, public  :: propagate_npt_mttk
  procedure, public  :: dump_baro_state

end type type_barostat

contains

! Initialise the thermostat
subroutine init_thermostat(th, inp, ndof, ke_atoms)

  ! passed variables
  class(type_thermostat), intent(inout)     :: th
  type(md_input), intent(in)                :: inp
  integer, intent(in)                       :: ndof
  real(double), intent(in)                  :: ke_atoms

  th%ndof = ndof
  th%ke_atoms = ke_atoms
  call th%get_temperature
  th%th_type = inp%thermo_type
  th%ensemble = inp%ensemble
  th%nat = inp%natoms
  th%dt = inp%dt
  th%T_ext = inp%T_ext
  th%iprint = inp%iprint
  th%tau_T = inp%tau_T
  th%k_B_md = one
  th%n_mts = inp%n_mts
  th%n_ys = inp%n_ys
  if (th%th_type == 'nhc') then
    th%n_nhc = inp%n_nhc
    call th%init_nhc(th%n_nhc, inp%nhc_mass, inp%cell_nhc_mass)
  end if
  call th%init_ys

end subroutine init_thermostat

! Initialise time steps for Yoshida-Suzuki integration
subroutine init_ys(th)

  ! passed variables
  class(type_thermostat), intent(inout)     :: th

  ! local variables
  integer                                   :: i, j, k, l, m
  real(double), dimension(:,:), allocatable :: psuz
  real(double)                              :: xnt

  allocate(psuz(th%n_ys,5))
  allocate(th%dt_ys(th%n_ys))

  do i=2,th%n_ys
    xnt = one/(two*real(i,double)-1)
    do j=1,5
      if (mod(j,3) == 0) then
        psuz(i,j) = one - four/(four-four**xnt)
      else
        psuz(i,j) = one/(four-four**xnt)
      end if
    end do
  end do

  ! Yoshida-Suzuki time steps
  k = 0
  select case(th%n_ys) ! The YS order
  case(1)
    th%dt_ys(1) = one
  case(3)
    th%dt_ys(1) = one/(two - two**(one/three))
    th%dt_ys(2) = one - two*th%dt_ys(1)
    th%dt_ys(3) = th%dt_ys(1)
  case(5)
    do i=1,5
      k = k+1
      th%dt_ys(k) = psuz(2,i)
    end do
  case(7)
    th%dt_ys(1) = -1.17767998417887_double
    th%dt_ys(2) = 0.235573213359357_double
    th%dt_ys(3) = 0.784513610477560_double
    th%dt_ys(4) = one - two*(th%dt_ys(1)+th%dt_ys(2)+th%dt_ys(3))
    th%dt_ys(5) = th%dt_ys(3)
    th%dt_ys(6) = th%dt_ys(2)
    th%dt_ys(7) = th%dt_ys(1)
  case(15)
    th%dt_ys(1) = 0.914844246229740_double
    th%dt_ys(2) = 0.253693336566229_double
    th%dt_ys(3) = -1.44485223686048_double
    th%dt_ys(4) = -0.158240635368243_double
    th%dt_ys(5) = 1.93813913762276_double
    th%dt_ys(6) = -1.96061023297549_double
    th%dt_ys(7) = 0.102799849391985_double
    th%dt_ys(8) = one - two*(th%dt_ys(1)+th%dt_ys(2)+th%dt_ys(3)+&
                             th%dt_ys(4)+th%dt_ys(5)+th%dt_ys(6)+&
                             th%dt_ys(7))
    th%dt_ys(9) = th%dt_ys(7)
    th%dt_ys(10) = th%dt_ys(6)
    th%dt_ys(11) = th%dt_ys(5)
    th%dt_ys(12) = th%dt_ys(4)
    th%dt_ys(13) = th%dt_ys(3)
    th%dt_ys(14) = th%dt_ys(2)
    th%dt_ys(15) = th%dt_ys(1)
  case(25)
    do j=1,5
      do i=1,5
        k = k+1
        th%dt_ys(k) = psuz(2,i)*psuz(3,j)
      end do
    end do
  case(125)
    do l=1,5
      do j=1,5
        do i=1,5
          k = k+1
          th%dt_ys(k) = psuz(2,i)*psuz(3,j)*psuz(4,l)
        end do
      end do
    end do
  case(625)
    do m=1,5
      do l=1,5
        do j=1,5
          do i=1,5
            k = k+1
            th%dt_ys(k) = psuz(2,i)*psuz(3,j)*psuz(4,l)*psuz(5,m)
          end do
        end do
      end do
    end do
  case default
    stop "Invalid Yoshida-Suzuki order"
  end select
  th%dt_ys = th%dt*th%dt_ys/th%n_mts
  deallocate(psuz)

end subroutine init_ys

subroutine get_temperature(th)

  ! passed variables
  class(type_thermostat), intent(inout)       :: th

  th%T_int = two*th%ke_atoms/th%k_B_md/real(th%ndof,double)
    write(*,'(6x,a,f16.8)') "Temperature = ", th%T_int

end subroutine get_temperature

! Isokinetic thermostat: maintain temperature by rescaling velocities by
! sqrt(T_int/T_ext) every tau_T steps
subroutine velocity_rescale(th, T_int, v)

  ! passed variables
  class(type_thermostat), intent(inout)       :: th
  real(double), intent(in)                    :: T_int
  real(double), dimension(:,:), intent(inout) :: v

  ! local variables
  real(double)  :: lambda

  lambda = sqrt(th%T_ext/T_int)
  v = v*lambda
  if (th%iprint > 1) write(*,'(4x,"Velocity rescale, lambda = ",f8.4)') lambda
end subroutine velocity_rescale

subroutine get_berendsen_thermo_sf(th, T_int)

  ! passed variables
  class(type_thermostat), intent(inout)       :: th
  real(double), intent(in)                    :: T_int

  th%lambda = sqrt(one + (th%dt/th%tau_T)*(th%T_ext/T_int - one))
  if (th%iprint > 1) then
    write(*,'(4x,"Weak coupling, thermostat lambda = ",f8.4)') th%lambda
  end if

end subroutine get_berendsen_thermo_sf

! Berendsen weak coupling thermostat
subroutine berendsen_thermo_propagate(th, v)

  ! passed variables
  class(type_thermostat), intent(inout)       :: th
  real(double), dimension(:,:), intent(inout) :: v

  v = v*th%lambda

end subroutine berendsen_thermo_propagate

! Initialise the NHC with n_nhc heat baths
subroutine init_nhc(th, n_nhc, m_nhc, m_nhc_cell)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: n_nhc
  real(double), dimension(:), intent(in)  :: m_nhc, m_nhc_cell

  th%n_nhc = n_nhc
  allocate(th%eta(n_nhc))
  allocate(th%v_eta(n_nhc))
  allocate(th%G_nhc(n_nhc))
  allocate(th%m_nhc(n_nhc))
  th%m_nhc = m_nhc
  if (th%cell_nhc) then
    allocate(th%eta_cell(n_nhc))
    allocate(th%v_eta_cell(n_nhc))
    allocate(th%G_nhc_cell(n_nhc))
    allocate(th%m_nhc_cell(n_nhc))
    th%m_nhc_cell = m_nhc_cell
  end if

  th%eta = zero
  if (th%n_nhc > 0) then
    th%m_nhc = one
    th%v_eta = sqrt(two*th%T_ext/th%m_nhc(1))
  else
    th%v_eta = zero
  end if
  th%G_nhc = zero
  call th%get_nhc_ke

end subroutine init_nhc

! get the force on thermostat k
subroutine update_G_nhc(th, k, box_ke)

  ! passed variables
  class(type_thermostat), intent(inout)   :: th
  integer, intent(in)                     :: k
  real(double), intent(in)                :: box_ke ! for P-T coupling

  if (th%cell_nhc) then
    if (k==1) then
      th%G_nhc(k) = 2*th%ke_atoms - th%ndof*th%k_B_md*th%T_ext
      th%G_nhc_cell(k) = 2*box_ke - th%T_ext*th%k_B_md*th%T_ext
    else
      th%G_nhc(k) = th%m_nhc(k-1)*th%v_eta(k-1)**2 - th%k_B_md*th%T_ext
      th%G_nhc_cell(k) = th%m_nhc_cell(k-1)*th%v_eta_cell(k-1)**2 - &
                         th%k_B_md*th%T_ext
    end if
    th%G_nhc(k) = th%G_nhc(k)/th%m_nhc(k)
    th%G_nhc_cell(k) = th%G_nhc_cell(k)/th%m_nhc_cell(k)
  else
    if (k == 1) then
      th%G_nhc(k) = 2*th%ke_atoms - th%ndof*th%k_B_md*th%T_ext + box_ke
    else
      th%G_nhc(k) = (th%m_nhc(k-1)*th%v_eta(k-1)**2 - th%k_B_md*th%T_ext)
    end if
    th%G_nhc(k) = th%G_nhc(k)/th%m_nhc(k)
  end if

  write(*,'(6x,"NHC: updating force, k = ",i2," G_nhc = ",f12.8)') k, th%G_nhc(k)

end subroutine update_G_nhc

! propagate eta (this is the same for all k)
subroutine propagate_eta(th, k, dt, dtfac)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer                               :: k
  real(double)                          :: dt, dtfac

  th%eta(k) = th%eta(k) + dtfac*dt*th%v_eta(k)
  if (th%cell_nhc) then
    th%eta_cell(k) = th%eta_cell(k) + dtfac*dt*th%v_eta(k)
  end if
  if (th%iprint > 1) then
    write(*,'(6x,"NHC: propagating eta,      k = ",i2," eta = ",f12.8)') k, th%eta(k)
    if (th%cell_nhc) then
      write(*,'(6x,"NHC: propagating eta_cell, k = ",i2," eta = ",f12.8)') k, th%eta_cell(k)
    end if
  end if

end subroutine propagate_eta

! propagate v_eta_k, simple shift step
subroutine propagate_v_eta_lin(th, k, dt, dtfac)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: k
  real(double), intent(in)              :: dt, dtfac

  th%v_eta(k) = th%v_eta(k) + dtfac*dt*th%G_nhc(k)
  if (th%cell_nhc) then
    th%v_eta_cell(k) = th%v_eta_cell(k) + dtfac*dt*th%G_nhc_cell(k)
  end if

  if (th%iprint > 1) then
    write(*,'(6x,"NHC: propagating v_eta linear,      k = ",i2," v_eta = ",f12.8)') k, th%v_eta(k)
    if (th%cell_nhc) then
      write(*,'(6x,"NHC: propagating v_eta_cell linear, k = ",i2," v_eta = ",f12.8)') k, th%v_eta_cell(k)
    end if
  end if

end subroutine propagate_v_eta_lin

! propagate v_eta_k, exponential shift step
subroutine propagate_v_eta_exp(th, k, dt, dtfac)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: k
  real(double), intent(in)              :: dt, dtfac

  th%v_eta(k) = th%v_eta(k)*exp(-dtfac*dt*th%v_eta(k+1))
  if (th%cell_nhc) then
  th%v_eta_cell(k) = th%v_eta_cell(k)*exp(-dtfac*dt*th%v_eta_cell(k+1))
  end if

  if (th%iprint > 1) then
    write(*,'(6x,"NHC: propagating v_eta exp,         k = ",i2," v_eta = ",f12.8)') k, th%v_eta(k)
    if (th%cell_nhc) then
      write(*,'(6x,"NHC: propagating v_eta_cell exp,    k = ",i2," v_eta = ",f12.8)') k, th%v_eta_cell(k)
    end if
  end if

end subroutine propagate_v_eta_exp

! Propagate the NHC
subroutine propagate_nvt_nhc(th, dt, v)

  ! passed variables
  class(type_thermostat), intent(inout)       :: th
  real(double), intent(in)                    :: dt
  real(double), dimension(:,:), intent(inout) :: v

  ! local variables
  integer       :: i, j, k
  real(double)  :: F_sys
  real(double)  :: dtys ! Yoshida-Suzuki time step
  real(double)  :: v_sfac ! velocity scaling factor
  real(double)  :: expfac

  if (th%iprint > 1) write(*,'(4x,a)') "Thermostat: propagating NHC dt/2 update"

  v_sfac = one
  call th%update_G_nhc(1, zero)
  do i=1,th%n_mts ! MTS loop
    do j=1,th%n_ys ! Yoshida-Suzuki loop
      ! Reverse part of expansion: update forces and thermostat velocities
      do k=th%n_nhc,1,-1
        if (k==th%n_nhc) then
          ! call th%update_G_nhc(k, zero)
          call th%propagate_v_eta_lin(k, dt, quarter)
        else
          ! Trotter expansion to avoid sinh singularity (see MTTK paper)
          call th%propagate_v_eta_exp(k, dt, eighth)
          call th%propagate_v_eta_lin(k, dt, quarter)
          ! call th%update_G_nhc(k, zero)
          call th%propagate_v_eta_exp(k, dt, eighth)
        end if
      end do
      expfac = exp(-half*dt*th%v_eta(1))
      v_sfac = v_sfac*expfac
      th%ke_atoms = th%ke_atoms*expfac**2
      call th%update_G_nhc(1, zero)
      ! update the nhc "positions" eta
      do k=1,th%n_nhc
        call th%propagate_eta(k, dt, half)
      end do
      ! Forward part of expansion: update forces and thermostat velocities
      call th%update_G_nhc(1, zero)
      do k=1,th%n_nhc
        if (k<th%n_nhc) then
          ! Trotter expansion to avoid sinh singularity (see MTTK paper)
          call th%propagate_v_eta_exp(k, dt, eighth)
          if (k /= 1) call th%update_G_nhc(k, zero)
          call th%propagate_v_eta_lin(k, dt, quarter)
          call th%propagate_v_eta_exp(k, dt, eighth)
        else
          call th%propagate_v_eta_lin(k, dt, quarter)
        end if
      end do
    end do ! Yoshida-Suzuki loop
  end do ! MTS loop

  v = v_sfac*v

end subroutine propagate_nvt_nhc

! Compute the "NHC energy" for monitoring the conserved quantity (KE + PE + NHC energy)
subroutine get_nhc_ke(th)

  ! passed variables
  class(type_thermostat), intent(inout) :: th

  if (th%n_nhc > 0) then
    th%e_nhc = half*sum(th%m_nhc*th%v_eta**2)
    th%e_nhc = th%e_nhc + th%ndof*th%k_B_md*th%T_ext*th%eta(1)
    th%e_nhc = th%e_nhc + th%k_B_md*th%T_ext*sum(th%eta(2:))
  else
    th%e_nhc = zero
  end if

  if (th%iprint > 1) then
    write(*,'(6x,a,e16.8)') "Updating nhc energy: e_nhc = ", th%e_nhc
  end if
  
end subroutine get_nhc_ke

! Debugging routine: dump the state of the thermostat
subroutine dump_thermo_state(th, step, funit)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: step, funit

  ! local variables
  integer                               :: i
  character(40)                         :: fmt


  write(fmt,'("(a8,",i4,"e16.4)")') th%n_nhc
  write(funit,'("step   ",i12)') step
  if (th%th_type == 'nhc') then
    write(funit,fmt) "eta:    ", th%eta
    write(funit,fmt) "v_eta:  ", th%v_eta
    write(funit,fmt) "G_nhc:  ", th%G_nhc
    write(funit,'("e_nhc:  ",e16.4)') th%e_nhc
  else if (th%th_type == 'berendsen') then
    write(funit,'("lambda: ",e16.4)') th%lambda
  end if
  write(funit,*)

end subroutine dump_thermo_state

subroutine init_barostat(baro, inp, ndof, ke_atoms, init_cell)

  use cell, only: type_cell

  ! Passed variables
  class(type_barostat), intent(inout)       :: baro
  type(md_input), intent(in)                :: inp
  type(type_cell), intent(in)               :: init_cell
  integer, intent(in)                       :: ndof
  real(double), intent(in)                  :: ke_atoms

  ! Local variables
  integer       :: i
  real(double)  :: m_box

  baro%ndof = ndof
  baro%ke_atoms = ke_atoms
  baro%iprint = inp%iprint
  baro%baro_type = inp%baro_type
  baro%ensemble = inp%ensemble
  baro%dt = inp%dt
  baro%nat = inp%natoms
  baro%P_ext = inp%P_ext
  baro%tau_P = inp%tau_P
  baro%V = init_cell%V
  baro%V_ref = baro%V
  baro%h = init_cell%h
  baro%h_0 = baro%h
  baro%bulkmod = inp%bulkmod

  baro%k_B_md = one
  baro%odnf = one + three/baro%ndof

  ! initialise MTTK
  baro%v_h = zero
  baro%m_box = inp%box_mass
  baro%ke_box = zero
  baro%ident = zero
  do i=1,3
    baro%ident(i,i) = one
  end do
  baro%stress_ext  = baro%ident*baro%P_ext
  baro%onfm = (one/baro%ndof)*baro%ident

end subroutine init_barostat

subroutine get_berendsen_baro_sf(baro, dt)

  ! Passed variables
  class(type_barostat), intent(inout)       :: baro
  real(double), intent(in)                  :: dt

  if (baro%iprint > 1) then
    write(*,'(6x,a)') "Getting berendsen barostat scaling factor"
  end if
  select case(baro%baro_type)
  case('berendsen')
    baro%mu = (one - (dt/baro%tau_P)*(baro%P_ext-baro%P_int)/baro%bulkmod)**third
    baro%mu_h = baro%ident*baro%mu
    if (baro%iprint > 1) then
      write(*,'(6x,a,f12.6)') "mu   = ", baro%mu
    end if
  case('ortho-berendsen')
    baro%mu_h = baro%ident - (dt/baro%tau_P)*(baro%stress_ext-baro%stress)/baro%bulkmod
    if (baro%iprint > 1) then
      write(*,'(6x,a,3f12.6)') "mu_h = ", baro%mu_h(1,:)
      write(*,'(6x,a,3f12.6)') "       ", baro%mu_h(2,:)
      write(*,'(6x,a,3f12.6)') "       ", baro%mu_h(3,:)
    end if
  end select
end subroutine get_berendsen_baro_sf

subroutine berendsen_baro_propagate(baro, r, h)

  ! Passed variables
  class(type_barostat), intent(inout)         :: baro
  real(double), dimension(:,:), intent(inout) :: r
  real(double), dimension(3,3), intent(inout) :: h

  ! local variables
  integer                                     :: i

  if (baro%iprint > 1) write(*,'(6x,a)') "Berendsen: propagating box and ions"

  select case(baro%baro_type)
  case('berendsen')
    baro%h = h*baro%mu
    h = baro%h
  case('ortho-berendsen')
  case('full-berendsen')
    baro%h = matmul(h, baro%mu_h)
    h = baro%h
  end select

  do i=1,baro%nat
    select case(baro%baro_type)
    case('berendsen')
      r(i,:) = r(i,:)*baro%mu
    case('ortho-berendsen')
    case('full-berendsen')
    end select
  end do

end subroutine berendsen_baro_propagate

! Update the box forces
subroutine update_G_box(baro, total_ke)

  ! Passed variables
  class(type_barostat), intent(inout)       :: baro
  real(double), intent(in)                  :: total_ke

  if (baro%iprint > 1) write(*,'(6x,a)') "MTTK: updating box force"

  select case(baro%baro_type)
  case('iso-mttk')
    baro%G_eps = (baro%odnf*two*total_ke + &
                 three*(baro%P_int - baro%P_ext)*baro%V)/baro%m_box
    if (baro%iprint > 1) then
      write(*,'(6x,a,f12.6)') "G_eps = ", baro%G_eps
    end if
  case('ortho-mttk')
  case('mttk')
    baro%G_h = (baro%ident*total_ke/baro%ndof + &
                baro%V*(baro%stress - baro%ident*baro%stress_ext))/baro%m_box
    if (baro%iprint > 1) then
      write(*,'(6x,a,3f12.6)') "G_h   = ", baro%G_h(1,:)
      write(*,'(6x,a,3f12.6)') "      = ", baro%G_h(2,:)
      write(*,'(6x,a,3f12.6)') "      = ", baro%G_h(3,:)
    end if
  end select

end subroutine update_G_box

subroutine propagate_eps_lin(baro, dt, dtfac)

  ! Passed variables
  class(type_barostat), intent(inout)   :: baro
  real(double), intent(in)              :: dt
  real(double), intent(in)              :: dtfac

  baro%eps = baro%eps + dtfac*dt*baro%v_eps
  if (baro%iprint > 1) write(*,'(6x,a,e16.8)') "MTTK: propagating eps linear shift: eps = ", baro%eps
  
end subroutine propagate_eps_lin

subroutine propagate_v_eps_lin(baro, dt, dtfac)

  ! Passed variables
  class(type_barostat), intent(inout)   :: baro
  real(double), intent(in)              :: dt
  real(double), intent(in)              :: dtfac

  baro%v_eps = baro%v_eps + dtfac*dt*baro%G_eps
  if (baro%iprint > 1) write(*,'(6x,a,f12.8)') "MTTK: propagating v_eps linear shift: v_eps = ", baro%v_eps
  
end subroutine propagate_v_eps_lin

subroutine propagate_v_eps_exp(baro, dt, dtfac, v_eta_1)

  ! Passed variables
  class(type_barostat), intent(inout)   :: baro
  real(double), intent(in)              :: dt
  real(double), intent(in)              :: dtfac
  real(double), intent(in)              :: v_eta_1

  baro%v_eps = baro%v_eps*exp(-dtfac*dt*v_eta_1)
  if (baro%iprint > 1) write(*,'(6x,a,f12.8)') "MTTK: propagating v_eps exp factor: v_eps = ", baro%v_eps
  
end subroutine propagate_v_eps_exp

subroutine propagate_r_ions(baro, dt, dtfac, v, r)

  ! Passed variables
  class(type_barostat), intent(inout)         :: baro
  real(double), intent(in)                    :: dt, dtfac
  real(double), dimension(:,:), intent(in)    :: v
  real(double), dimension(:,:), intent(inout) :: r

  ! Local variables
  integer                                     :: i, j
  real(double)                                :: sinhx_x, rscale, vscale, &
                                                 exp_v_eps, eps
  real(double), dimension(3)                  :: vp
  real(double), dimension(3,3)                :: rfac_m, vfac_m

  select case(baro%baro_type)
  case('iso-mttk')
    exp_v_eps = exp(dt*dtfac*baro%v_eps)
    sinhx_x = poly_sinhx_x(dtfac*dt*baro%v_eps)
    rscale = exp_v_eps**2
    vscale = exp_v_eps*sinhx_x*dt
    r = r*rscale + v*vscale
  case('ortho-mttk')
  case('mttk')
    rfac_m = matmul(baro%c_g, baro%I_e)
    vfac_m = matmul(baro%c_g, baro%I_s)
    do i=1,baro%nat
      vp = matmul(v(i,:), vfac_m)*dt
      r(i,:) = matmul(r(i,:), rfac_m) + vp
      r(i,:) = matmul(r(i,:), transpose(baro%c_g))
    end do
  end select

  if (baro%iprint > 1) &
    write(*,'(4x,a)') "MTTK: propagating particle positions"

end subroutine propagate_r_ions

subroutine propagate_v_ions(baro, dt, dtfac, v_eta_1, v)

  ! Passed variables
  class(type_barostat), intent(inout)         :: baro
  real(double), intent(in)                    :: dt, dtfac, v_eta_1
  real(double), dimension(:,:), intent(inout) :: v

  ! local variables
  integer                                     :: i
  real(double)                                :: vfac
  real(double), dimension(3,3)                :: vfac_m

  select case(baro%baro_type)
  case('iso-mttk')
    vfac = exp(-dtfac*dt*(v_eta_1 + baro%odnf*baro%v_eps))
    v = v*vfac
    baro%ke_atoms = baro%ke_atoms*vfac**2
  case('ortho-mttk')
  case('mttk')
    vfac_m = matmul(baro%c_g, baro%I_e)
    vfac_m = matmul(vfac_m, transpose(baro%c_g))
    do i=1,baro%nat
      v(i,:)  = matmul(v(i,:), vfac_m)
    end do
  end select

  if (baro%iprint > 1) &
    write(*,'(6x,a)') "MTTK: propagating particle velocities"

end subroutine propagate_v_ions

subroutine propagate_box(baro, h, V)

  ! Passed variables
  class(type_barostat), intent(inout)         :: baro
  real(double), dimension(3,3), intent(inout) :: h
  real(double), intent(inout)                 :: V

  ! Local variables
  integer                                     :: i
  real(double)                                :: V_new, V_old, hscale
  real(double), dimension(3,3)                :: htemp

  select case(baro%baro_type)
  case('iso-mttk')
    V_old = baro%V
    V_new = exp(three*baro%eps)*baro%V_ref
    hscale = (V_new/V_old)**third
    baro%h_scale = baro%ident*hscale
    baro%h = baro%h*hscale
    baro%V = V_new
    V = baro%V
    h = baro%h
  case('ortho-mttk')
  case('mttk')
    baro%h = transpose(baro%h)
    htemp = matmul(baro%h, baro%c_g)
    baro%h = matmul(htemp, baro%I_e)
    htemp = matmul(baro%h, transpose(baro%c_g))
    baro%h = htemp
    V = baro%V
    h = baro%h
  end select

  if (baro%iprint > 1) write(*,'(6x,a)') "MTTK: propagating box"

end subroutine propagate_box

! Get the eigenvalues/vectors of v_h(0) + Tr[v_h(0)]/N_f + v_eta(1)
subroutine diag_vbox(baro, v_eta_1)

  ! Passed variables
  class(type_barostat), intent(inout)      :: baro
  real(double), intent(in)                  :: v_eta_1

  ! Local variables
  real(double), dimension(:), allocatable   :: work
  real(double)                              :: tr_vh, sinhx_x
  integer                                   :: i, lwork, info
  lwork = 27
  allocate(work(lwork))

  baro%c_g = zero
  tr_vh = zero
  do i=1,3
    tr_vh = tr_vh + baro%v_h(i,i)
  end do
  baro%c_g = baro%v_h + (tr_vh/baro%ndof + v_eta_1)*baro%ident
  ! diagonalise
  call dsyev('V', 'U', 3, baro%c_g, 3, baro%lambda, work, lwork, info)

  if (baro%iprint > 1) write(*,'(6x,a)') "MTTK: Computing eigenvalues of v_h"

end subroutine diag_vbox

subroutine get_Ie(baro, dt, dtfac)

  ! Passed variables
  class(type_barostat), intent(inout)       :: baro
  real(double), intent(in)                  :: dt, dtfac

  ! local variables
  integer                                   :: i

  baro%I_e = zero
  do i=1,3
    baro%I_e(i,i) = exp(baro%lambda(i)*dt*dtfac)
  end do

end subroutine get_Ie

subroutine get_Is(baro, dt, dtfac)

  ! Passed variables
  class(type_barostat), intent(inout)       :: baro
  real(double), intent(in)                  :: dt, dtfac

  ! local variables
  integer                                   :: i
  real(double)                              :: sinhx_x

  baro%I_s = zero
  do i=1,3
    sinhx_x = poly_sinhx_x(baro%lambda(i)*dt*dtfac)
    baro%I_s(i,i) = exp(baro%lambda(i)*dt*dtfac)*sinhx_x
  end do

end subroutine get_Is

! Get the box kinetic energy
subroutine get_box_ke(baro)

  ! Passed variables
  class(type_barostat), intent(inout)      :: baro

  ! local variables
  real(double), dimension(3,3)              :: temp
  integer                                   :: i

  select case(baro%baro_type)
  case('iso-mttk')
    baro%ke_box = half*baro%m_box*baro%v_eps**2
  case('ortho-mttk')
    baro%ke_box = half*baro%m_box*sum(baro%v_Q**2)
  case('mttk')
    temp = matmul(transpose(baro%v_h), baro%v_h)
    baro%ke_box = zero
    do i=1,3
      baro%ke_box = baro%ke_box + temp(i,i)
    end do
    baro%ke_box = baro%m_box*baro%ke_box
  end select

end subroutine get_box_ke

subroutine update_stress(baro, static_stress, kinetic_stress)

  ! passed variables
  class(type_barostat), intent(inout)       :: baro
  real(double), dimension(3,3), intent(in)  :: static_stress, kinetic_stress

  baro%static_stress = static_stress
  baro%kinetic_stress = kinetic_stress
end subroutine update_stress

! Compute the stress tensor from the virial tensor from get_force_and_energy
! and the kinetic energy from get_kinetic_energy
subroutine get_stress_and_pressure(baro)

  ! passed variables
  class(type_barostat), intent(inout)       :: baro

  ! local variables
  integer                         :: i

  baro%stress = -(one/baro%V)*(baro%kinetic_stress + baro%static_stress)
  if (baro%iprint > 1) then
    write(*,*)
    write(*,'(4x,a36,4x,a36)') "kinetic stress", "static stress"
    do i=1,3
      write(*,'(4x,3f12.4,4x,3f12.4)') baro%kinetic_stress(i,:), &
                                       baro%static_stress(i,:)
    end do
    write(*,*)
  end if

  baro%P_int = zero
  do i=1,3
    baro%P_int = baro%P_int + baro%stress(i,i)
  end do
  baro%P_int = baro%P_int*third

  write(*,'("  Pressure         = ",f16.8)') baro%P_int
  write(*,'(2x,a)') "Stress tensor:"
  do i=1,3
    write(*,'(4x,3f16.8)') baro%stress(i,:)
  end do

end subroutine get_stress_and_pressure

! linear shift in box velocity
subroutine propagate_v_h_1(baro, dt, dtfac)

  ! Passed variables
  class(type_barostat), intent(inout)      :: baro
  real(double), intent(in)                 :: dt, dtfac

  if (baro%iprint > 1) write(*,'(6x,a)') "MTTK: propagating v_h linear shift"
  baro%v_h = baro%v_h + dtfac*dt*baro%G_h

end subroutine propagate_v_h_1

! exponential shift in box velocity from thermostat coupling
subroutine propagate_v_h_2(baro, dt, dtfac, v_eta)

  ! Passed variables
  class(type_barostat), intent(inout)      :: baro
  real(double), intent(in)                  :: dt, dtfac, v_eta

  if (baro%iprint > 1) write(*,'(6x,a)') "MTTK: propagating v_h exp factor"
  baro%v_h = baro%v_h*exp(-dtfac*dt*v_eta)

end subroutine propagate_v_h_2

! Update the box
subroutine propagate_h(baro, dt, dtfac, h)

  ! Passed variables
  class(type_barostat), intent(inout)         :: baro
  real(double), intent(in)                    :: dt, dtfac
  real(double), dimension(3,3), intent(inout) :: h

  ! local variables
  integer                                     :: i
  real(double), dimension(3,3)                :: htemp

  if (baro%iprint > 1) write(*,'(4x,a)') "MTTK: propagating h"

  baro%h = transpose(h)
  ! htemp = matmul(baro%h, baro%c_g)
  ! baro%h = matmul(htemp, baro%I_e)
  ! htemp = matmul(baro%h, transpose(baro%c_g))

  htemp = mmult3x3_loop(baro%c_g, baro%h)
  baro%h = mmult3x3_loop(htemp, baro%I_e)
  htemp = mmult3x3_loop(transpose(baro%c_g), baro%h)
  ! htemp = mmult3x3(baro%h, baro%c_g)
  ! h = mmult3x3(htemp, baro%I_e)
  ! htemp = mmult3x3(h, transpose(baro%c_g))
  ! baro%h = htemp
  h = htemp
  baro%h = htemp

end subroutine propagate_h

! Integrate the thermostat and barostat
subroutine propagate_npt_mttk(baro, th, dt, ke_atoms, v)

  ! Passed variables
  class(type_barostat), intent(inout)         :: baro
  class(type_thermostat), intent(inout)       :: th
  real(double), intent(in)                    :: dt, ke_atoms
  real(double), dimension(:,:), intent(inout) :: v

  ! Local variables
  integer                                    :: i_ys, i_mts, i_nhc, i_atom
  real(double)                               :: v_eta_couple, v_sfac

  if (baro%iprint > 1) write(*,'(4x,a)') "MTTK: integrating thermo/barostat"

  v_sfac = one
  call baro%get_box_ke
  call th%update_G_nhc(1, baro%ke_box)
  call baro%update_G_box(ke_atoms)

  do i_mts=1,th%n_mts ! Multiple time step loop
    do i_ys=1,th%n_ys ! Yoshida-Suzuki loop
      ! Update thermostat velocities
      call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
      do i_nhc=th%n_nhc-1,-1
        call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), eighth)
        call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
        call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), eighth)
      end do

      ! Update box velocities
      if (th%cell_nhc) then
        v_eta_couple = th%v_eta_cell(1)
      else
        v_eta_couple = th%v_eta(1)
      end if
      select case(baro%baro_type)
      case('iso-mttk')
        call baro%propagate_v_eps_exp(th%dt_ys(i_ys), eighth, v_eta_couple)
        call baro%propagate_v_eps_lin(th%dt_ys(i_ys), quarter)
        call baro%propagate_v_eps_exp(th%dt_ys(i_ys), eighth, v_eta_couple)
      case('ortho-mttk')
      case('mttk')
      end select

      ! propagate atomic velocities
      call baro%propagate_v_ions(th%dt_ys(i_ys), half, v_eta_couple, v)
      th%ke_atoms = baro%ke_atoms

      if (baro%baro_type == 'iso-mttk' .or. baro%baro_type == 'ortho-mttk' &
          .or. baro%baro_type == 'mttk') then
        call baro%update_G_box(ke_atoms)
      end if

      ! Propagate thermostat positions
      do i_nhc=1,th%n_nhc
        call th%propagate_eta(i_nhc, th%dt_ys(i_ys), half)
      end do

      ! Update box velocities
      if (th%cell_nhc) then
        v_eta_couple = th%v_eta_cell(1)
      else
        v_eta_couple = th%v_eta(1)
      end if
      select case(baro%baro_type)
      case('iso-mttk')
        call baro%propagate_v_eps_exp(th%dt_ys(i_ys), eighth, v_eta_couple)
        call baro%propagate_v_eps_lin(th%dt_ys(i_ys), quarter)
        call baro%propagate_v_eps_exp(th%dt_ys(i_ys), eighth, v_eta_couple)
      case('ortho-mttk')
      case('mttk')
      end select
      call baro%propagate_v_ions(th%dt_ys(i_ys), half, th%v_eta(1), v)

      call baro%get_box_ke
      call th%update_G_nhc(1, baro%ke_box)
      ! Update thermostat velocities
      do i_nhc=1,th%n_nhc-1
        call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), eighth)
        if (i_nhc /= 1) call th%update_G_nhc(i_nhc, baro%ke_box)
        call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
        call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), eighth)
      end do
      call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
    end do ! Yoshida-Suzuki loop
  end do ! Multiple time step loop


end subroutine propagate_npt_mttk

! Debugging routine: dump the state of the barostat
subroutine dump_baro_state(baro, step, funit)

  ! Passed variables
  class(type_barostat), intent(inout)        :: baro
  integer, intent(in)                        :: step, funit

  ! local variables
  integer                                    :: i

  write(funit,'("step   ",i12)') step
  write(funit,'("box ke ",e12.4)') baro%ke_box
  write(funit,'(8x,a)') "Lattice vectors"
  do i=1,3
    write(funit,'(4x,3f12.4)') baro%h(i,:)
  end do
  write(funit,'(8x,a)') "Lattice scaling factors"
  do i=1,3
    write(funit,'(4x,3f12.4)') baro%h_scale(i,:)
  end do
  write(funit,'(8x,a,29x,a)') "Virial KE", "Stress"
  do i=1,3
    write(funit,'(4x,3e12.4,4x,3e12.4)') baro%kinetic_stress(i,:), &
                                         baro%stress(i,:)
  end do
  select case(baro%baro_type)
  case('iso-mttk')
    write(funit,'("eps:   ",e16.4)') baro%eps
    write(funit,'("v_eps: ",e16.4)') baro%v_eps
    write(funit,'("ke_box:",e16.4)') baro%ke_box
    write(funit,'("G_eps: ",e16.4)') baro%G_eps
  case('ortho-mttk')
  case('mttk')
    write(funit,'(8x,a)') 'v_h'
    do i=1,3
      write(funit,'(4x,3e12.4)') baro%v_h(i,:)
    end do
    write(funit,'(8x,a,12x,a)') 'lambda', 'c_g'
    do i=1,3
      write(funit,'(4x,e12.4,6x,3e12.4)') baro%lambda(i), baro%c_g(i,:)
    end do
    write(funit,'(8x,a)') "G_h"
    do i=1,3
      write(funit,'(4x,3e12.4)') baro%G_h(i,:)
    end do
  end select
    write(funit,'("P_int: ",e16.4)') baro%P_int
    write(funit,'("V:     ",e16.4)') baro%V
  write(funit,*)

end subroutine dump_baro_state

end module md_control