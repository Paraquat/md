program md

use cmdline
use md_module

implicit none

type(type_cmdline)                  :: opts

call opts%get_args()


end program md
