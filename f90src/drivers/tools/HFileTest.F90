program MAIN
  use HistFileMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  character(len=*), parameter :: prognm=__FILE__
  real(r8), pointer :: ptr_1d(:)

  call hist_addfld1d(fname='TSOI',units='oC',avgflag='A',long_name='soil temperature',ptr_col=ptr_1d)
  call hist_addfld1d(fname='OMC', units='gC m^-^3',avgflag='A',long_name='microbial carbon',ptr_col=ptr_1d)
  call hist_addfld1d(fname='OMN', units='gN m^-^3',avgflag='A',long_name='microbial nitrogen',ptr_col=ptr_1d)
  call hist_addfld1d(fname='OMP', units='gP m^-^3',avgflag='A',long_name='microbial phosphorus',ptr_col=ptr_1d)

  call hist_printflds()
end program main
