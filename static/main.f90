program md

use cmdline
use datatypes
use constants
use cell
use static_module
use input_module

implicit none

type(type_cmdline)                  :: opts
type(static_input)                  :: inp
type(type_cell)                     :: init_cell
type(type_static)                   :: static
type(type_pairpotential)            :: pp

integer   :: i
character(40) :: inpfile

inpfile = 'md.in'

call opts%get_args()

! read parameters from md.in
call inp%read_input(inpfile)
call init_cell%read_cell(inp%sfile)

call pp%init_pp(opts%ppfile, opts%shift)

call static%init_static(inp, init_cell, pp)
call static%run_static

end program md
