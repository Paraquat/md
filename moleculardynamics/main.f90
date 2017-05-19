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

call mdr%init_md(init_cell, pp, pp%ns, opts%shift)
call mdr%get_force_and_energy
write(*,*) "Energy = ", mdr%pe_t
write(*,*) "Force = ", mdr%f_t(1,1)

end program md
