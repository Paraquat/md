module progress
! Progress bar for a loop; writes to stderr
! initialise with init_progress: requires total number of loops (maxloop) and 
! text to display above progress bar
! update with update_progress in loop: requires current loop number

use datatypes

implicit none

private
public :: type_progress

type type_progress
  integer               :: maxloops
  character(18)         :: barstr
  character(10)         :: bar

  contains
    procedure :: init_progress
    procedure :: update_progress
end type type_progress

contains

subroutine init_progress(p, maxloops, displaytext)

  class(type_progress)  :: p

  character(80)         :: displaytext
  integer               :: maxloops

  p%maxloops = maxloops
  p%bar = "          "

  write(0,*) displaytext

end subroutine init_progress

subroutine update_progress(p, loop)

  class(type_progress)  :: p

  integer, intent(in)   :: loop
  integer               :: i
  real(double)          :: progress, pc

  progress = real(loop,double)/real(p%maxloops,double)
  do i=1,nint(progress*10)
    p%bar(i:i) = "#"
  end do

  write(p%barstr,"(a,i3,a3,a10,a)") '\r',nint(progress*100),"% [",p%bar,"]"
  if (loop .ge. p%maxloops) then
    write(0,'(a)',advance='yes') p%barstr
  else
    write(0,'(a)',advance='no') p%barstr
  end if
end subroutine update_progress

end module progress
