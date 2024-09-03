program ClimReader
  use ClimReadMod
  use abortutils, only : endrun
  USE EcoSIMCtrlMod, ONLY : LVERB
  use EcoSIMCtrlMod, only : clm_hour_file_in
implicit none

  character(len=40) :: progname
  integer :: num_args
  integer :: iyear
  type(atm_forc_type)  :: atmf
  integer :: L,I
  integer :: irec
  character(len=20) :: buf

  CALL GETARG(0,progname)

  num_args = iargc()

  if (num_args<2)then
    call endrun(msg='Please use code as: '//trim(progname)//' infile iyear')
  endif
  call GETARG(1,clm_hour_file_in)
  call GETARG(2,buf);read(buf,*)iyear

  lverb=.true.
  irec=1
  L=2
  call ReadClimNC(iyear,irec,L,atmf)

  print*,atmf
end program ClimReader
