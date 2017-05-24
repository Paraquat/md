program md

use cmdline
use datatypes
use constants
use cell
use static_module

implicit none

type(type_cmdline)                  :: opts
type(type_cell)                     :: init_cell
type(type_static)                   :: static
type(type_pp)                       :: pp

integer   :: i

call opts%get_args()

call pp%init_pp_from_file(opts%ppfile, opts%shift)
call init_cell%read_cell(opts%sfile)

call static%init_static(init_cell, pp, opts%cart, opts%shift)
call static%get_force_and_energy

end program md
