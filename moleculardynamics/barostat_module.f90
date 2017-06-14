module barostat_module

use datatypes
use constants
use maths, only: exp_sinhx_x

implicit none

type type_barostat

  character(40)                 :: baro_type
  real(double)                  :: P_ext  ! applied pressure
  real(double)                  :: W_h    ! box mass
  real(double), dimension(3,3)  :: v_h    ! box velocity
  real(double), dimension(3,3)  :: stress ! stress tensor
  real(double), dimension(3,3)  :: ke_stress ! ke contribution to stress tensor
  real(double), dimension(3,3)  :: press  ! full pressure tensor including ke
  real(double), dimension(3,3)  :: G_h    ! box force
  real(double), dimension(3,3)  :: ident  ! 3x3 identity matrix
  real(double), dimension(3,3)  :: c_g    ! eigenvectors 
  real(double), dimension(3)    :: lambda ! eigenvalues
  integer                       :: iprint, nat, ndof
  real(double)                  :: ke_box
  real(double)                  :: k_B_md

contains

  procedure :: init_barostat
  procedure :: init_xs
  procedure :: get_lambda_cg
  procedure :: get_lambda_sfac
  procedure :: get_box_ke
  procedure :: update_G_h
  procedure :: propagate_h
  procedure :: propagate_v_h_1
  procedure :: propagate_v_h_2
  procedure :: propagate_v_sys
  procedure :: propagate_r_sys
  procedure :: vVerlet_r_h_npt

end type type_barostat

contains

  subroutine init_barostat(baro, baro_type, P_ext, nat, ndof, box_mass, iprint)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    character(40), intent(in)                 :: baro_type
    real(double), intent(in)                  :: P_ext, box_mass
    integer, intent(in)                       :: nat, ndof, iprint

    ! Local variables
    integer       :: i
    real(double)  :: W_h

    baro%baro_type = baro_type
    baro%nat = nat
    baro%ndof = ndof
    baro%iprint = iprint

    baro%v_h = zero
    baro%k_B_md = one

    baro%W_h = box_mass

    baro%ident = zero
    do i=1,3
      baro%ident(i,i) = one
    end do

    if (baro%baro_type == 'mttk') call baro%init_xs

  end subroutine init_barostat

  ! General subroutine for initialiseing extended-system/Lagrangian barostats
  ! (Andersen, Parrinello-Rahma, MTTK, ...)
  subroutine init_xs(baro)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro

  end subroutine init_xs

  ! Get the eigenvalues/vectors of v_h(0) + Tr[v_h(0)]/N_f + v_eta(1)
  subroutine get_lambda_cg(baro, v_eta)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: v_eta

    ! Local variables
    real(double), dimension(:), allocatable   :: work
    real(double)                              :: tr_vh
    integer                                   :: i, lwork, info

    lwork = 27
    allocate(work(lwork))

    baro%c_g = zero
    tr_vh = zero
    do i=1,3
      tr_vh = tr_vh + baro%v_h(i,i)
    end do

    ! baro%c_g = baro%v_h + (tr_vh/baro%ndof + v_eta)*baro%ident
    baro%c_g = baro%v_h + tr_vh/baro%ndof + v_eta
    ! diagonalise
    call dsyev('V', 'U', 3, baro%c_g, 3, baro%lambda, work, lwork, info)

  end subroutine get_lambda_cg

  subroutine get_lambda_sfac(baro, dt, dtfac, rscale, vscale)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: dt, dtfac
    real(double), dimension(3), intent(out)   :: rscale, vscale

    ! local variables
    integer                                   :: i
    real(double)                              :: exp_eig, sinhx_x

    do i=1,3
      exp_eig = exp(dtfac*dt*baro%lambda(i))
      sinhx_x = exp_sinhx_x(dtfac*dt*baro%lambda(i))
      rscale(i) = exp_eig**2
      vscale(i) = exp_eig*dt*sinhx_x
    end do

  end subroutine get_lambda_sfac

  ! Get the box kinetic energy
  subroutine get_box_ke(baro)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro

    ! local variables
    real(double), dimension(3,3)              :: temp
    integer                                   :: i

    temp = matmul(transpose(baro%v_h), baro%v_h)
    baro%ke_box = zero
    do i=1,3
      baro%ke_box = baro%ke_box + temp(i,i)
    end do
    baro%ke_box = baro%W_h*baro%ke_box

  end subroutine get_box_ke

  ! Update the box forces
  subroutine update_G_h(baro, total_ke, virial_ke, volume, T_ext, Q_1, G_nhc_1)

    ! Passed variables
    class(type_barostat), intent(inout)       :: baro
    real(double), dimension(3,3), intent(in)  :: virial_ke
    real(double), intent(in)                  :: total_ke
    real(double), intent(in)                  :: volume
    real(double), intent(in)                  :: T_ext
    real(double), intent(in)                  :: Q_1 ! mass of NHC k=1
    real(double), intent(out)                 :: G_nhc_1 ! force on NHC k=1

    baro%G_h = (baro%ident*total_ke/baro%ndof + virial_ke + &
                volume*(baro%press - baro%ident*baro%P_ext))/baro%W_h
    ! The force on the first NHC heat bath is required for thermostat coupling
    G_nhc_1 = (total_ke + baro%ke_box - baro%ndof*baro%k_B_md*T_ext)/Q_1

  end subroutine update_G_h

  ! Update the box
  subroutine propagate_h(baro, rscale, h)

    ! Passed variables
    class(type_barostat), intent(inout)        :: baro
    real(double), dimension(3), intent(in)      :: rscale
    real(double), dimension(3,3), intent(inout) :: h

    ! local variables
    integer                                   :: i, j
    real(double), dimension(3,3)              :: htemp

    htemp = matmul(transpose(baro%c_g), h)
    do i=1,3
      do j=1,3
        htemp(j,i) = htemp(j,i)*rscale(j)
      end do 
    end do
    h = matmul(baro%c_g, htemp)

  end subroutine propagate_h

  ! linear shift in box velocity
  subroutine propagate_v_h_1(baro, dt, dtfac)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: dt, dtfac

    baro%v_h = baro%v_h + dtfac*dt*baro%G_h

  end subroutine propagate_v_h_1

  ! exponential shift in box velocity from thermostat coupling
  subroutine propagate_v_h_2(baro, dt, dtfac, v_eta)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: dt, dtfac, v_eta

    baro%v_h = baro%v_h*exp(-dtfac*dt*v_eta)

  end subroutine propagate_v_h_2

  ! coupling box and particle velocities
  subroutine propagate_v_sys(baro, dt, dtfac, v_eta_1, v)

    ! Passed variables
    class(type_barostat), intent(inout)         :: baro
    real(double), intent(in)                    :: dt, dtfac
    real(double), intent(in)                    :: v_eta_1
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    real(double), dimension(3)                :: vscale, vtemp
    integer                                   :: i, j

    call baro%get_lambda_cg(v_eta_1)

    do i=1,3
      vscale(i) = exp(-baro%lambda(i)*dt*dtfac)
    end do

    do i=1,baro%nat
      do j=1,3
        vtemp(j) = dot_product(v(i,:), baro%c_g(:,j))
      end do
      vtemp = vtemp*vscale
      do j=1,3
        v(i,j) = dot_product(vtemp, baro%c_g(:,j))
      end do
    end do

  end subroutine propagate_v_sys

  ! coupling box and particle positions. This is a little complicated and
  ! contains a polynomial expansion of sinh(x)/x, so I will do it in one
  ! subroutine for now, as in the MTTK paper.
  subroutine propagate_r_sys(baro, dt, rscale, vscale, r,  v)

    ! Passed variables
    class(type_barostat), intent(inout)        :: baro
    real(double), intent(in)                    :: dt
    real(double), dimension(:,:), intent(inout) :: r
    real(double), dimension(:,:), intent(in)    :: v
    real(double), dimension(3), intent(in)      :: rscale, vscale

    ! local variables
    integer                                     :: i, j
    real(double), dimension(3)                  :: u, uv
    real(double)                                :: exp_eig, sinhx_x

    do i=1,baro%nat
      do j=1,3
        u(j) = dot_product(r(i,:), baro%c_g(:,j))
        uv(j) = dot_product(v(i,:), baro%c_g(:,j))
        u(j) = u(j)*rscale(j) + uv(j)*vscale(j)
      end do
      do j=1,3
        r(i,j) = dot_product(u, baro%c_g(j,:))
      end do
    end do

  end subroutine propagate_r_sys

  ! Propagate the ionic positions and box of the extended system. Temperature
  ! coupling determined by v_eta
  subroutine vVerlet_r_h_npt(baro, dt, h, v, r, v_eta)

    ! Passed variables
    class(type_barostat), intent(inout)        :: baro
    real(double), intent(in)                    :: dt
    real(double), intent(in)                    :: v_eta ! thermostat velocity
    real(double), dimension(3,3), intent(inout) :: h
    real(double), dimension(:,:), intent(inout) :: r
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    real(double), dimension(3)                  :: rscale, vscale

    ! Preceded by velocity Verlet dt/2 step
    call baro%get_lambda_cg(v_eta)
    call baro%get_lambda_sfac(dt, half, rscale, vscale)
    call baro%propagate_r_sys(dt, rscale, vscale, r, v)
    call baro%propagate_h(rscale, h)
    ! Proceded by velocity Verlet dt/2 step

  end subroutine vVerlet_r_h_npt

end module barostat_module
