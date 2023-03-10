program MAIN
  use HistFileMod
  implicit none
  character(len=*), parameter :: prognm=__FILE__


  call hist_addfld1d('TSOI','oC')
  call hist_addfld1d('OMC','gC m^-^3')
  call hist_addfld1d('OMN','gN m^-^3')
  call hist_addfld1d('OMP','gP m^-^3')

  call hist_printflds()
end program main
