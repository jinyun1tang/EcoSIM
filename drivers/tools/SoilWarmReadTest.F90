PROGRAM SoilWarmReadTest

  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod     , only : etimer,warming_exp
  use ecosim_time_mod   , only : getdow
  use abortutils        , only : endrun
  use ClimReadMod       , only : read_soil_warming_Tref
  use ncdio_pio
  use netcdf
  use PerturbationMod
  implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer :: year

  warming_exp='4K;1m;Blodget.ctrl.ecosim.h1.xxxx-01-00-00000.nc;2014/01/01:2018/12/31'  

  call config_soil_warming(warming_exp)

  do year=2013,2020
   if(check_Soil_Warming(year,1))then
!     call read_soil_warming_Tref(year)
   endif
  enddo
end PROGRAM SoilWarmReadTest
