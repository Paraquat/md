module trajectory
! container for a MD trajectory. Stores positions r, velocities v, forces f
! and lattice vectors lat

use cell
use progress
use datatypes

implicit none

type type_traj
  type(type_cell) :: initcell
  integer         :: nsteps, natoms, nspec, nt
  real(double)    :: T

  real(double), allocatable, dimension(:,:,:) :: r
  real(double), allocatable, dimension(:,:,:) :: v
  real(double), allocatable, dimension(:,:,:) :: f
  real(double), allocatable, dimension(:,:,:) :: h

  character(2), allocatable, dimension(:) :: species

  contains
    procedure :: init_traj
    procedure :: update
    procedure :: init_vasp
    procedure :: read_outcar
    procedure :: init_xsf
    procedure :: read_xsf
!    procedure :: read_md
end type type_traj

contains

! initialise
subroutine init_traj(traj, nsteps, natoms)

  class(type_traj), intent(inout) :: traj

  integer, intent(in)             :: nsteps, natoms

  traj%nt = 0
  traj%natoms = natoms
  traj%nsteps = nsteps+1

  allocate(traj%r(traj%nsteps,traj%natoms,3))
  allocate(traj%v(traj%nsteps,traj%natoms,3))
  allocate(traj%f(traj%nsteps,traj%natoms,3))
  allocate(traj%h(traj%nsteps,3,3))
end subroutine init_traj

! update trajectory with a new cell object
subroutine update(traj, cell)

  class(type_traj), intent(inout) :: traj

  class(type_cell), intent(in)    :: cell

  traj%r(traj%nt,:,:) = cell%r
  traj%nt = traj%nt + 1
end subroutine update

! zeroth MD step is the intial configuration in POSCAR
subroutine init_vasp(traj, infilename, nspec)

  class(type_traj), intent(inout) :: traj


  character(40), intent(in)       :: infilename
  integer, intent(in)             :: nspec

  call traj%initcell%read_poscar(infilename, nspec)

  traj%nspec = nspec
  traj%r(1,:,:) = traj%initcell%r
  traj%h(1,:,:) = traj%initcell%h
  traj%species = traj%initcell%species
  traj%nt = 0

end subroutine init_vasp

! Read a MD trajectory (r and f only) from OUTCAR file
subroutine read_outcar(traj)

  class(type_traj), intent(inout) :: traj

  type(type_progress)             :: pbar
  character(40)                   :: infilename
  character(100)                  :: buf, junk, displaytext
  integer                         :: fh, stat, i

  infilename = "OUTCAR"
  displaytext = "Reading trajectory from OUTCAR..."
  fh = 102

  open(fh, file=infilename, status='old')

  call pbar%init_progress(traj%nsteps,displaytext)
  outer: do
    if (traj%nt .gt. traj%nsteps) then
      write(*,*) "nt exceeded nsteps. Stopping."
      write(*,*) traj%nt, traj%nsteps
      exit outer
    end if
    read(fh,'(a)',iostat=stat) buf
    if (stat /= 0) exit outer
    if (index(buf,'direct lattice vectors') .gt. 0) then
      if (traj%nt .eq. 0) then ! the first instance is the POSCAR vectors
        traj%nt = traj%nt+1
        cycle outer
      end if
      traj%nt = traj%nt + 1
      call pbar%update_progress(traj%nt)
      inner1: do i=1,3
        read(fh,*,iostat=stat) traj%initcell%h(i,:)
        if (stat /= 0) exit outer
      end do inner1
      traj%h(traj%nt,:,:) = traj%initcell%h
      call traj%initcell%invert_lat()
    else if (index(buf,'POSITION') .gt. 0) then
      read(fh,'(a)',iostat=stat) junk
      if (stat /= 0) exit outer
      inner2: do i=1,traj%natoms
        read(fh,*,iostat=stat) traj%initcell%r(i,:), traj%f(traj%nt,i,:)
        if (stat /= 0) exit outer
        ! because for SOME reason the OUTCAR contains Cartesian coordinates
        traj%r(traj%nt,i,:) = traj%initcell%cart2frac(traj%initcell%r(i,:))
      end do inner2
    end if
  end do outer

  close(fh)
end subroutine read_outcar

subroutine init_xsf(traj, infilename, nspec)

  class(type_traj), intent(inout) :: traj


  character(40), intent(in)       :: infilename
  integer, intent(in)             :: nspec

  call traj%initcell%read_xsf(infilename)

  traj%nspec = nspec
  traj%r(1,:,:) = traj%initcell%r
  traj%h(1,:,:) = traj%initcell%h
  traj%species = traj%initcell%species
  traj%nt = 0
end subroutine init_xsf

! read a .xsf file
subroutine read_xsf(traj, infilename, nat, nspec, speclist)

  class(type_traj), intent(inout)     :: traj

  integer, intent(in)                 :: nspec, nat
  character(2), dimension(10), intent(in) :: speclist

  type(type_progress)                 :: pbar
  type(type_cell)                     :: cell

  character(40)                       :: infilename, junk
  integer                             :: i, j, nframe, step

  infilename = trim(infilename)
  allocate(cell%r(nat,3),cell%rcart(nat,3),cell%species(nat))
  cell%nspec = nspec
  cell%h = 0.

  open(101, file=infilename, status='old')
  read(101,*) junk, nframe
  read(101,*) junk
  do i=1,nframe
    read(101,*) junk, step
    read(101,*) cell%h(1,:)
    read(101,*) cell%h(2,:)
    read(101,*) cell%h(3,:)
    read(101,*) junk
    read(101,*) junk
    call cell%invert_lat()

    do j=1,nat
      read(101,*) cell%species(j), cell%rcart(j,:)
      traj%species(j) = cell%species(j)
      traj%r(i,j,:) = cell%cart2frac(cell%rcart(j,:))
      traj%nt = traj%nt + 1
    end do
  end do

  close(101)
  call cell%count_species()
end subroutine read_xsf
!subroutine read_md
!end subroutine read_md

end module trajectory
