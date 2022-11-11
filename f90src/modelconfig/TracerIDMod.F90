module TracerIDMod

  use MiniMathMod, only : addone
implicit none
  save
  CHARACTER(LEN=*), private, PARAMETER :: MOD_FILENAME=__FILE__
  integer, parameter :: idg_CO2=1
  integer, parameter :: idg_CH4=2
  integer, parameter :: idg_O2 =3
  integer, parameter :: idg_N2 =4
  integer, parameter :: idg_N2O=5
  integer, parameter :: idg_H2 =6
  integer, parameter :: idg_NH3=7
  integer, parameter :: idg_NH3B=8

  integer :: ids_NH4,ids_NH4B
  integer :: ids_NO3,ids_NO3B
  integer :: ids_NO2,ids_NO2B
  integer :: ids_H1PO4,ids_H1PO4B
  integer :: ids_H2PO4,ids_H2PO4B
  integer :: ids_DOC,ids_DON,ids_DOP,ids_ODA
  integer :: idg_beg,idg_end
  integer :: ids_beg,ids_end
  integer :: ids_nut_beg,ids_nuts_beg,ids_nuts_end
  contains

  subroutine InitTracerIDs(lsalt_model)
  implicit none
  logical, intent(in) :: lsalt_model

  idg_beg=1; ids_beg=1
! for better array manipulation of land-atmosphere exchange,
! banded NH3 is considered as (potential) gas too. However, there
! is no banded gas NH3 concentration.

  idg_end=idg_NH3B;

  ids_nuts_beg=idg_NH3;  !the first nutrient tracer, including band
  ids_end=idg_end   !initalize the solute counter
  ids_NH4B=addone(ids_end);ids_NO3B=addone(ids_end);
  ids_NO2B=addone(ids_end);ids_H1PO4B=addone(ids_end)
  ids_H2PO4B=addone(ids_end)

  ids_NH4=addone(ids_end);
  ids_NO3=addone(ids_end);
  ids_NO2=addone(ids_end);
  ids_H1PO4=addone(ids_end);
  ids_H2PO4=addone(ids_end);

  ids_nut_beg=ids_NH4;  !the first non-band non-gaseous nutrient tracer
  ids_nuts_end=ids_H2PO4;!the last non-band nutrient tracer
  if(lsalt_model)then

  else

  endif
  end subroutine InitTracerIDs
end module TracerIDMod
