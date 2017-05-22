program md

use cmdline
use datatypes
use constants
use cell
use md_module

implicit none

type(type_cmdline)                  :: opts
type(type_cell)                     :: init_cell
type(type_md)                       :: mdr
type(type_pp)                       :: pp

integer   :: i

call opts%get_args()

call pp%init_pp_from_file(opts%ppfile, opts%shift)
call init_cell%read_cell(opts%sfile)

call mdr%init_md(init_cell, pp, opts%cart, opts%ensemble, pp%ns, opts%dt, &
                 opts%T_ext, opts%shift, opts%comv)
call mdr%md_run(1,opts%nsteps)

end program md
