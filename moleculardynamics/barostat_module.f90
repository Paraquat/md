module barostat_module

use datatypes
use constants
use rng

implicit none

type type_barostat

contains

  procedure :: weak_coupling

end type type_barostat

contains

  subroutine weak_coupling
  end subroutine weak_coupling

end module barostat_module
