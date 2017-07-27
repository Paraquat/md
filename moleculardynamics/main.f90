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

! Read arguments from command line
call opts%get_args(args)
if (args .eqv. .true.) then
  nsteps = opts%nsteps
  call pp%init_pp(opts%ppfile, opts%shift)
  call init_cell%read_cell(opts%sfile)

  call mdr%init_md(init_cell, pp, opts%cart, opts%ensemble, pp%ns, opts%dt, &
                   opts%T_ext, opts%init_distr, opts%shift, opts%comv, &
                   opts%ttype, opts%tau_T, opts%n_nhc, opts%nhc_mass, &
                   opts%btype, opts%P_ext, opts%boxm, opts%tau_P)
  mdr%dump_freq = opts%dump_freq
  mdr%tau_T = opts%tau_T
else
  ! otherwise try to read parameters from md.in
  call inp%read_input(inpfile)
  nsteps = inp%nsteps
  call pp%init_pp(inp%ppfile, inp%shift)
  call init_cell%read_cell(inp%sfile)

  call mdr%init_md(init_cell, pp, opts%cart, inp%ensemble, pp%ns, inp%dt, &
                   inp%T_ext, inp%v_distr, inp%shift, inp%comv, &
                   inp%thermo_type, inp%tau_T, inp%n_nhc, inp%nhc_mass, &
                   inp%baro_type, inp%P_ext, inp%box_mass, inp%tau_P)
  mdr%dump_freq = inp%dump_freq
  mdr%tau_T = inp%tau_T
end if

call mdr%md_run(1,nsteps)

end program md
