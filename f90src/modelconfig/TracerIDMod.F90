module TracerIDMod

  use MiniMathMod, only : addone
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
  integer, parameter :: idg_NH3B=8

  integer :: ids_NH4,ids_NH4B
  integer :: ids_NO3,ids_NO3B
  integer :: ids_NO2,ids_NO2B
  integer :: ids_H1PO4,ids_H1PO4B
  integer :: ids_H2PO4,ids_H2PO4B
  integer :: ids_DOC,ids_DON,ids_DOP,ids_ODA
  integer :: idg_beg,idg_end
  integer :: ids_beg,ids_end
  contains

  subroutine InitTracerIDs(lsalt_model)
  implicit none
  logical, intent(in) :: lsalt_model

  idg_beg=1; ids_beg=1
  idg_end=idg_H2; ids_end=idg_end

  ids_NH4=addone(ids_end);ids_NH4B=addone(ids_end)
  ids_NO3=addone(ids_end);ids_NO3B=addone(ids_end)
  ids_NO2=addone(ids_end);ids_NO2B=addone(ids_end)
  ids_H1PO4=addone(ids_end);ids_H1PO4B=addone(ids_end)
  ids_H2PO4=addone(ids_end);ids_H2PO4B=addone(ids_end)

  if(lsalt_model)then

  else

  endif
  end subroutine InitTracerIDs
end module TracerIDMod
