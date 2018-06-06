program md

use cmdline
use input_module
use datatypes
use constants
use cell
use md_module

implicit none

type(type_cmdline)                  :: opts
type(md_input)                      :: inp
type(type_cell)                     :: init_cell
type(type_md)                       :: mdr
type(type_pairpotential)            :: pp
logical                             :: args
character(40)                       :: inpfile

integer       :: nsteps

inpfile = 'md.in'

! Initialisation
! Read arguments from command line
call opts%get_args(args)
! read parameters from md.in
call inp%read_input(inpfile)
nsteps = inp%nsteps
call pp%init_pp(inp%ppfile, inp%shift)
call init_cell%read_cell(inp%sfile)
init_cell%pbc_method = 'frac'

call mdr%init_md(inp, init_cell, pp)
mdr%dump_freq = inp%dump_freq
mdr%cp_freq = inp%cp_freq

! Run the main MD loop
call mdr%md_run(mdr%step,nsteps)

end program md
