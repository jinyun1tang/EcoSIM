module TracerIDMod

implicit none
  save
  CHARACTER(LEN=*), private, PARAMETER :: MOD_FILENAME=__FILE__
  integer, parameter :: idg_CO2=1
  integer, parameter :: idg_CH4=2
  integer, parameter :: idg_O2 =3
  integer, parameter :: idg_NH3=4
  integer, parameter :: idg_N2 =5
  integer, parameter :: idg_N2O=6
  integer, parameter :: idg_H2 =7

  integer :: idg_beg,idg_end

  contains

  subroutine InitTracerIDs(lsalt_model)
  implicit none
  logical, intent(in) :: lsalt_model

  idg_beg=1
  idg_end=idg_H2

  if(lsalt_model)then

  else

  endif
  end subroutine InitTracerIDs
end module TracerIDMod
