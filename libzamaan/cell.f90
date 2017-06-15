module cell

use vector
use datatypes
use constants

implicit none

private
public :: type_cell

type type_cell
  integer         :: nat
  integer         :: nspec

  real(double), allocatable, dimension(:,:) :: r, rcart, dt
  character(2), allocatable, dimension(:)   :: spec, species
  integer, allocatable, dimension(:)        :: spec_count, spec_int
  real(double), allocatable, dimension(:)   :: mass
  real(double), dimension(3,3)              :: h, h_inv
  real(double), dimension(6)                :: param

  contains
    procedure :: count_species
    procedure :: int_label_species
    procedure :: recell
    procedure :: uncell
    procedure :: volume
    procedure :: frac2cart
    procedure :: disp_frac2cart
    procedure :: disp_frac2cart_noshift
    procedure :: cell_frac2cart
    procedure :: invert_lat
    procedure :: cart2frac
    procedure :: mic
    procedure :: wrap_positions_cart
    procedure :: init
    procedure :: get_dt
    procedure :: read_xyz
    procedure :: read_shx
    procedure :: read_cell
    procedure :: read_poscar
    procedure :: read_lammps
    procedure :: read_xsf
    procedure :: write_poscar
    procedure :: write_cell
    procedure :: write_lammps
    procedure :: write_xsf
    procedure :: supercell
    procedure :: cut_ortho
end type type_cell

contains

subroutine count_species(p)

  class(type_cell), intent(inout)   :: p

  integer, parameter                :: max_spec=5
  integer                           :: i, j, nspecies
  integer, dimension(max_spec)      :: scount_temp
  character(2), dimension(max_spec) :: s_temp
  logical                           :: notfound

  s_temp = "z"
  scount_temp = 0
  nspecies = 0

  do i=1,p%nat
    notfound = .true.
    do j=1,max_spec
      if (p%species(i) .eq. s_temp(j)) then
        scount_temp(j) = scount_temp(j) + 1
        notfound = .false.
      end if
    end do
    if (notfound .eqv. .true.) then
      do j=1,max_spec
        if (s_temp(j) .eq. "z") then
          s_temp(j) = p%species(i)
          scount_temp(j) = 1
          nspecies = nspecies+1
          exit
        end if
      end do
    end if
  end do

  p%nspec = nspecies
  allocate(p%spec(p%nspec),p%spec_count(p%nspec))
  do i=1,p%nspec
    p%spec(i)=s_temp(i)
    p%spec_count(i) = scount_temp(i)
  end do

end subroutine count_species

! Give species integer label
subroutine int_label_species(p)

  class(type_cell), intent(inout)   :: p

  integer                           :: i, j

  allocate(p%spec_int(p%nat))

  do i=1,p%nat
    do j=1,p%nspec
      if (p%species(i) == p%spec(j)) then
        p%spec_int(i) = j
      end if
    end do
  end do

end subroutine int_label_species

! calculate lattice vectors from cell parameters
subroutine recell(p)

  class(type_cell), intent(inout)   :: p

  real(double)                      :: a, b, c, alp_deg, bet_deg, gam_deg, &
                                       cosa, cosb, cosg, sing

  a = p%param(1)
  b = p%param(2)
  c = p%param(3)

  alp_deg = p%param(4)
  bet_deg = p%param(5)
  gam_deg = p%param(6)

  if (p%param(4)-90.0_double < small) then
    cosa=0.0
  else
    cosa = cos(alp_deg*deg2rad)
  endif
  if (p%param(5)-90.0_double < small) then
    cosb = 0.0
  else
    cosb = cos(bet_deg*deg2rad)
  endif
  if (p%param(6)-90.0_double < small) then
    cosg = 0.0
    sing = 1.0
  else
    cosg = cos(gam_deg*deg2rad)
    sing = sin(gam_deg*deg2rad)
  endif

  p%h = 0.0
  p%h(1,1) = a
  p%h(2,1) = b*cosg
  p%h(2,2) = b*sing
  p%h(3,1) = c*cosb
  p%h(3,2) = c*(cosa - cosb*cosg)/sing
  p%h(3,3) = sqrt(c**2 - p%h(3,1)**2 - p%h(3,2)**2)
end subroutine recell

! calculate lattice parameters from cell vectors
subroutine uncell(p)

  class(type_cell), intent(inout)   :: p

  real(double)                      :: u, v, w

  p%param(1) = sqrt(sum(p%h(1,:)**2))
  p%param(2) = sqrt(sum(p%h(2,:)**2))
  p%param(3) = sqrt(sum(p%h(3,:)**2))

  u = dot_product(p%h(2,:), p%h(3,:))
  v = dot_product(p%h(1,:), p%h(3,:))
  w = dot_product(p%h(1,:), p%h(2,:))

  p%param(4) = rad2deg*acos(u/(p%param(2)*p%param(3)))
  p%param(5) = rad2deg*acos(v/(p%param(1)*p%param(3)))
  p%param(6) = rad2deg*acos(w/(p%param(1)*p%param(2)))
end subroutine uncell

! scalar triple product volume of cell
function volume(p) result(vol)

  class(type_cell), intent(inout)   :: p

  real(double), dimension(3)        :: cp
  real(double)                      :: vol

  cp = cross_product(p%h(2,:), p%h(3,:))
  vol = dot_product(p%h(1,:), cp)
end function volume

function frac2cart(p, v) result (r)

  class(type_cell), intent(inout)   :: p

  real(double), dimension(3), intent(in)    :: v
  real(double), dimension(3)                :: r, w

  integer                           :: i

  w=v
  do i=1,3
    if ( w(i) .gt. one ) then
      w(i) = w(i)-one
    else if ( w(i) .lt. zero) then
      w(i) = w(i) + one
    end if
  end do

  r = matmul(w,p%h)
end function frac2cart

! convert displacement from franctional to cartesian, shift coordinates
! outside the unit cell
function disp_frac2cart(p, v) result (r)

  class(type_cell), intent(inout)   :: p

  real(double), dimension(3), intent(in)    :: v
  real(double), dimension(3)                :: r, w

  integer                           :: i

  w=v
  do i=1,3
    if ( w(i) .gt. 0.5 ) then
      w(i) = w(i)-1.
    else if ( w(i) .lt. -0.5) then
      w(i) = w(i) + 1.
    end if
  end do

  r = matmul(w,p%h)
end function disp_frac2cart

! convert displacement from franctional to cartesian, with no shift for
! coordinates outside the unit cell
function disp_frac2cart_noshift(p, v) result (r)

  class(type_cell), intent(inout)   :: p

  real(double), dimension(3), intent(in)    :: v
  real(double), dimension(3)                :: r

  r = matmul(v,p%h)
end function disp_frac2cart_noshift

! Convert the cell from fractional to Cartesian
subroutine cell_frac2cart(p)

  class(type_cell), intent(inout)   :: p

  integer :: i

  if (allocated(p%rcart) .eqv. .false.) allocate(p%rcart(p%nat,3))

  do i=1,p%nat
    p%rcart(i,:) = p%frac2cart(p%r(i,:))
  end do

end subroutine cell_frac2cart

! Invert basis matrix
subroutine invert_lat(p)

  class(type_cell), intent(inout)   :: p

  real(double), dimension(3,3)   :: dum
  real(double), dimension(3)     :: work
  integer, dimension(3)     :: ipiv
  integer                   :: i, info

  ipiv = 0
  i = 3
  dum = p%h
  call DGETRF(i,i,dum,i,ipiv,info)
  call DGETRI(i,dum,i,ipiv,work,i,info)
  p%h_inv = dum

end subroutine invert_lat

! convert cartesian coordinates to fractional
function cart2frac(p,v) result(w)

  class(type_cell), intent(inout)   :: p

  real(double), dimension(3), intent(in)    :: v
  real(double), dimension(3)                :: w
  integer                           :: i

  w=matmul(v,p%h_inv)
  do
    do i=1,3
      if (w(i) .ge. 1.) then
        w(i)=w(i)-1.
      else if (w(i) .lt. 0.) then
        w(i)=w(i)+1.
      end if

      if (abs(w(i)-1.) .lt. 1.E-10) then
        w(i) = 0.
      end if
    end do
    if (maxval(w) .lt. 1.01 .and. minval(w) .gt. -1.E-12) then
      exit
    end if
  end do

end function cart2frac

! Minimum image convention periodic boundary conditions
! (only works for orthorhombic cells)
function mic(p, coord1, coord2) result(relative)

  class(type_cell), intent(inout)         :: p
  real(double), dimension(3), intent(in)  :: coord1, coord2
  real(double), dimension(3)              :: relative

  integer                            :: i

  relative = coord2 - coord1
  do i=1,3
    relative(i) = relative(i) - p%h(i,i)*nint(relative(i)/p%h(i,i))
  end do
end function mic

! Wrap atoms into the unit cell (orthrhombic cell, Cartesian coordinates only)
subroutine wrap_positions_cart(p)

  class(type_cell), intent(inout)   :: p

  integer                           :: i, j

  do i=1,p%nat
    do j=1,3
      if (p%rcart(i,j) > p%h(j,j)) then
        p%rcart(i,j) = p%rcart(i,j) - p%h(j,j)
      else if (p%rcart(i,j) < zero) then
        p%rcart(i,j) = p%rcart(i,j) + p%h(j,j)
      end if
    end do
  end do

end subroutine wrap_positions_cart

! Initialise cell from data in memory
subroutine init(p, nat, nspec, h, r, species)

  class(type_cell), intent(inout)   :: p

  real(double), dimension(:,:), intent(in)        :: r
  real(double), dimension(3,3), intent(in)        :: h
  character(2), dimension(:), intent(in)  :: species
  integer, intent(in)                     :: nat, nspec

  p%nat = nat
  p%nspec = nspec
  allocate(p%r(nat,3),p%rcart(nat,3),p%species(nat))
  p%r = r
  p%h = h
  p%species = species

  call p%count_species
  call p%int_label_species
  call p%cell_frac2cart

end subroutine init

! get the distance table
subroutine get_dt(p)

  class(type_cell), intent(inout)     :: p

  integer                                 :: iat, jat
  real(double), dimension(3)              :: r_ij, r_ij_cart
  real(double)                            :: d

  p%dt=0.

  do iat=1,p%nat
    do jat=iat+1,p%nat
      r_ij=p%r(jat,:)-p%r(iat,:)
      r_ij_cart=p%disp_frac2cart(r_ij)
      d=sqrt(sum(r_ij_cart**2))
      p%dt(iat,jat)=d
      p%dt(jat,iat)=d
    end do
  end do

end subroutine get_dt

! read xyz file (cell parameters in the comment line)
subroutine read_xyz(p, infilename)

  class(type_cell), intent(inout)   :: p

  character(40)                     :: infilename
  integer                           :: i

  infilename = trim(infilename)

  open(101, file=infilename, status='old')
  read(101,*) p%nat
  read(101,*) p%param

  allocate(p%rcart(p%nat,3), p%species(p%nat))
  do i = 1, p%nat
    read(101,*) p%species(i), p%rcart(i,:)
  end do

  close(101)
end subroutine read_xyz

! read shelx .res file (from AIRSS output)
subroutine read_shx(p, infilename)

  class(type_cell), intent(inout)   :: p

  character(40)                     :: infilename, id, junk, sg
  integer                           :: q1, q2, q3, i
  real(double)                      :: volume, energy, q4

  infilename = trim(infilename)

  open(101, file=infilename, status='old')
  read(101,*) junk, id, q1, volume, energy, q2, q3, p%nat, sg
  allocate(p%r(p%nat,3), p%species(p%nat))
  read(101,*) junk, q4, p%param
  read(101,*) junk, q1
  read(101,*) junk
  do i = 1, p%nat
    read(101,*) p%species(i), q1, p%r(i,:), q4
  end do
  close(101)
end subroutine read_shx

! read castep .cell file
subroutine read_cell(p, infilename)

  class(type_cell), intent(inout)   :: p

  character(40)                       :: infilename, s
  character(80)                       :: line
  integer                             :: i, j, reason
  real(double)                        :: m

  infilename = trim(infilename)

  p%nat = 0
  open(101, file=infilename, status='old')
  ! read lattice vectors
  outer1: do
    read(101,'(a)',iostat=reason) line
    if (reason > 0) then
      write(*,*) "Read error - this shouldn't happen!"
    else if (reason < 0) then
      exit outer1
    else
      if (line(1:19) .eq. "%BLOCK lattice_cart") then
        inner1: do i = 1,3
          read(101,*) p%h(i,:)
        end do inner1
        rewind(unit=101)
        exit outer1
      end if
    end if
  end do outer1

  rewind(unit=101)
  ! read ionic positions
  outer2: do
    read(101,'(a)',iostat=reason) line
    if (reason > 0) then
      write(*,*) "Read error - this shouldn't happen!"
    else if (reason < 0) then
      exit outer2
    else
      if (line(1:21) .eq. "%BLOCK positions_frac") then
        if (p%nat .eq. 0) then
          inner2: do
            read(101,'(a)') line
            if (line(1:24) .eq. "%ENDBLOCK positions_frac") then
              rewind(unit=101)
              exit inner2
            else
              p%nat = p%nat + 1
            end if
          end do inner2
        else
          allocate(p%species(p%nat), p%r(p%nat,3))
          do i = 1, p%nat
            read(101,*) p%species(i), p%r(i,:)
          end do
        end if
      end if
    end if
  end do outer2


  call p%cell_frac2cart
  call p%count_species
  call p%int_label_species
  allocate(p%mass(p%nspec))
  rewind(unit=101)

  ! read masses
  outer3: do
    read(101,'(a)',iostat=reason) line
    if (reason > 0) then
      write(*,*) "Read error - this shouldn't happen!"
    else if (reason < 0) then
      exit outer3
    else
      if (trim(line(1:19)) .eq. "%BLOCK species_mass") then
        inner3: do i=1,p%nspec
          read(101,*) s, m
          do j=1,p%nspec
            if (trim(s) == p%spec(j)) p%mass(j) = m
          end do
        end do inner3
        rewind(unit=101)
        exit outer3
      end if
    end if
  end do outer3

  close(101)

end subroutine read_cell

! read vasp poscar file
subroutine read_poscar(p, infilename, nspec)

  class(type_cell), intent(inout)   :: p

  character(40)                       :: infilename, junk
  integer                             :: i, j, acount, nspec
  real(double)                        :: sfac

  infilename = trim(infilename)
  p%nspec = nspec
  allocate(p%spec(nspec), p%spec_count(nspec))

  open(101, file=infilename, status='old')
  read(101,*) junk
  read(101,*) sfac
  read(101,*) p%h(1,:)
  read(101,*) p%h(2,:)
  read(101,*) p%h(3,:)
  read(101,*) p%spec(:)
  read(101,*) p%spec_count(:)

  p%h = p%h*sfac
  p%nat = sum(p%spec_count)
  allocate(p%species(p%nat), p%r(p%nat,3))
  read(101,*) junk
  acount = 1
  do i = 1,p%nspec
    do j = 1,p%spec_count(i)
      read(101,*) p%r(acount,:)
      p%species(acount) = p%spec(i)
      acount = acount + 1
    end do
  end do
  close(101)
  call p%invert_lat()
end subroutine read_poscar

! read lammps dump file
! Species MUST be sorted by atom id via script line:
! dump_modify dumpID sort id
subroutine read_lammps(p, infilename, nspec, speclist)

  class(type_cell), intent(inout)     :: p
  integer, intent(in)                 :: nspec
  character(2), dimension(10), intent(in) :: speclist

  character(40)                       :: infilename, junk
  integer                             :: i, junkint, ispec
  real(double)                        :: timestep, junkreal

  infilename = trim(infilename)
  p%nspec = nspec
  p%h = 0.

  open(101, file=infilename, status='old')
  read(101,*) junk
  read(101,*) timestep
  read(101,*) junk
  read(101,*) p%nat
  read(101,*) junk
  read(101,*) junkreal, p%h(1,1)
  read(101,*) junkreal, p%h(2,2)
  read(101,*) junkreal, p%h(3,3)
  read(101,*) junk
  call p%invert_lat()

  allocate(p%species(p%nat), p%r(p%nat,3), p%rcart(p%nat,3))

  do i=1,p%nat
    read(101,*) junkint, ispec, p%rcart(i,:)
    p%r(i,:) = p%cart2frac(p%rcart(i,:))
    p%species(i) = speclist(ispec)
  end do

  close(101)
  call p%count_species()
end subroutine read_lammps

! read .xsf file
subroutine read_xsf(p, infilename)

  class(type_cell), intent(inout)   :: p

  character(40)                       :: junk, infilename, testfile
  integer                             :: i, nspec, junkint

  testfile = 'test.res'
  infilename = trim(infilename)
  p%nspec = nspec

  open(101, file=infilename, status='old')
  read(101,*) junk
  read(101,*) junk
  read(101,*) p%h(1,:)
  read(101,*) p%h(2,:)
  read(101,*) p%h(3,:)
  read(101,*) junk
  read(101,*) p%nat, junkint

  allocate(p%species(p%nat), p%r(p%nat,3), p%rcart(p%nat,3))
  call p%invert_lat()

  do i=1,p%nat
    read(101,*) p%species(i), p%rcart(i,:)
    p%r(i,:) = p%cart2frac(p%rcart(i,:))
  end do

  close(101)

  call p%count_species()
end subroutine read_xsf

! write cell to poscar
subroutine write_poscar(p, filename)

  class(type_cell), intent(inout)   :: p

  character(40)                       :: filename
  integer                             :: i

  open(101, file=filename, status='replace')
  write(101,*) filename
  write(101,*) 1.0
  write(101,*) p%h(1,:)
  write(101,*) p%h(2,:)
  write(101,*) p%h(3,:)
  write(101,*) p%spec
  write(101,*) p%spec_count
  write(101,*) "Direct"
  do i = 1, p%nat
    write(101,*) p%r(i,:)
  end do
  close(101)
end subroutine write_poscar

! write cell to CASTEP cell file
subroutine write_cell(p, filename)

  class(type_cell), intent(inout)   :: p

  character(40)                       :: filename
  integer                             :: i

  open(101, file=filename, status='replace')
  write(101,*) "%BLOCK lattice_abc"
  write(101,*) p%h(1,:)
  write(101,*) p%h(2,:)
  write(101,*) p%h(3,:)
  write(101,*) "%ENDBLOCK lattice_abc"
  write(101,*)
  write(101,*) "%BLOCK positions_frac"
  do i = 1, p%nat
    write(101,*) p%species(i), p%r(i,:)
  end do
  write(101,*) "%ENDBLOCK positions_frac"
  close(101)
end subroutine write_cell

! write cell to LAMMPS input file
subroutine write_lammps(p, filename, comment)

  class(type_cell), intent(inout)     :: p
  character(80), intent(in)           :: comment

  character(40)                       :: filename
  integer                             :: i, j, s

  open(101, file=filename, status='replace')
  write(101,'(a)') trim(comment)
  write(101,'(i8,a6)') p%nat,' atoms'
  write(101,'(i8,a11)') p%nspec, ' atom types'
  write(101,*)
  write(101,'(2f8.4,a8)') 0.0, p%h(1,1), ' xlo xhi'
  write(101,'(2f8.4,a8)') 0.0, p%h(2,2), ' ylo yhi'
  write(101,'(2f8.4,a8)') 0.0, p%h(3,3), ' zlo zhi'
  write(101,*)
  write(101,'(a)') 'Masses'
  write(101,*)
  do i=1,p%nspec
    write(101,'(i4)') i
  end do
  write(101,*)
  write(101,'(a)') 'Atoms'
  write(101,*)
  do i=1,p%nat
    do j=1,p%nspec
      if (p%species(i) .eq. p%spec(j)) then
        s = j
        exit
      end if
    end do
    write(101,'(i8,i4,f6.3,3f14.8)') i, s, 0.0, p%disp_frac2cart(p%r(i,:))
  end do
  close(101)
end subroutine write_lammps

! Write structure to .xsf file (for MD visualisation)
subroutine write_xsf(p, iunit, md, step, nsteps)

  ! passed variales
  class(type_cell), intent(inout)   :: p
  logical, intent(in)               :: md
  integer, intent(in)               :: iunit, step, nsteps

  ! local variables
  integer                           :: i

  if (step == 1) write(iunit,'("ANIMSTEP ",i8)') nsteps
  write(iunit,'(a)') "CRYSTAL"
  write(iunit,'("PRIMVEC   ", i8)') step
  do i=1,3
    write(iunit,'(3f14.8)') p%h(i,:)
  end do
  write(iunit,'("PRIMCOORD ", i8)') step
  write(iunit,'(2i8)') p%nat, 1
  do i=1,p%nat
    write(iunit,'(a4,3f16.8)') p%species(i), p%rcart(i,:)
  end do

end subroutine write_xsf

! Construct nx x ny x nz supercell
subroutine supercell(p, sc_dim, psuper)

  class(type_cell), intent(inout)   :: p
  integer, dimension(3), intent(in) :: sc_dim ! (/nx, ny, nz/)

  class(type_cell), intent(out)     :: psuper

  integer                           :: natsuper, i, j, k, m, n, ncell
  real(double), dimension(3)                :: trans
  real(double), dimension(3,3)              :: latsuper
  real(double), dimension(:,:), allocatable :: rsuper
  character(2), dimension(:), allocatable   :: specsuper

  natsuper=p%nat*sc_dim(1)*sc_dim(2)*sc_dim(3)

  latsuper(1,:) = p%h(1,:)*sc_dim(1)
  latsuper(2,:) = p%h(2,:)*sc_dim(2)
  latsuper(3,:) = p%h(3,:)*sc_dim(3)

  allocate(rsuper(natsuper,3),specsuper(natsuper))

  ncell = 0
  do i=1, sc_dim(1)
    do j=1, sc_dim(2)
      do k=1, sc_dim(3)
        trans=(/i,j,k/)
        do m=1,p%nat
          n=ncell*p%nat+m
          rsuper(n,:)=(p%r(m,:)+trans)/sc_dim
          specsuper(n)=p%species(m)
        end do
        ncell=ncell+1
      end do
    end do
  end do

  call psuper%init(natsuper, p%nspec, latsuper, rsuper, specsuper)
  deallocate(rsuper,specsuper)

end subroutine supercell

! calculate the distance cutoff for an orthogoanl unit cell ONLY
subroutine cut_ortho(p, cutoff)

  class(type_cell), intent(inout)   :: p

  real(double), intent(out)                 :: cutoff

  cutoff=min(p%h(1,1), p%h(2,2), p%h(3,3))/2
end subroutine cut_ortho

end module cell
