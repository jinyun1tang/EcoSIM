program managementReader
  use ReadManagementMod
  use EcoSIMCtrlMod, only : soil_mgmt_in, lverb
  use abortutils, only : endrun
  USE FertilizerDataType
  use FlagDataType
implicit none

  character(len=40) :: progname
  integer :: num_args
  integer :: iyear

  CALL GETARG(0,progname)

  num_args = iargc()

  if (num_args<1)then
    call endrun(msg='Please use code as: '//trim(progname)//' infile')
  endif
  call GETARG(1,soil_mgmt_in)

  CALL InitFertilizerData

  call InitFlagData
  iyear=1
  lverb=.true.

  call ReadManagementFiles(iyear)

end program managementReader
