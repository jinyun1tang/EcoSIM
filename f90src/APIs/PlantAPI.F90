module PlantAPI

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use ExtractMod   , only : extract
  use grosubMod    , only : grosub
  use HfuncMod     , only : hfunc
  use UptakeMod    , only : uptake
  use timings      , only : start_timer, end_timer
  use EcoSIMCtrlDataType
implicit none

  private
  public :: PlantModel
  contains

  subroutine PlantModel(I,J,NHW,NHE,NVN,NVS)


  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: t1

333   FORMAT(A8)

!   UPDATE PLANT PHENOLOGY IN 'HFUNC'
!
  if(lverb)WRITE(*,333)'HFUNC'
  call start_timer(t1)
  CALL HFUNC(I,J,NHW,NHE,NVN,NVS)
  call end_timer('HFUNC',t1)
!
!   CALCULATE CANOPY CO2 UPTAKE AT FULL TURGOR, CANOPY WATER POTENTIAL,
!   HYDRAULIC AND STOMATAL RESISTANCES,AND CANOPY ENERGY BALANCE IN 'UPTAKE'
!   CALCULATE ROOT UPTAKE OF WATER, OXYGEN, NH4, NO3 AND PO4 IN 'UPTAKE'
!
  if(lverb)WRITE(*,333)'UPTK'
!    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL UPTAKE(I,J,NHW,NHE,NVN,NVS)
  call end_timer('UPTK',t1)
!
!   CALCULATE CANOPY CO2 UPTAKE AT AMBIENT TURGOR, AUTOTROPHIC AND GROWTH
!   RESPIRATION, PLANT C ALLOCATION, CANOPY AND ROOT GROWTH IN 'GROSUB'
!
  if(lverb)WRITE(*,333)'GRO'
!    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL GROSUB(I,J,NHW,NHE,NVN,NVS)
  call end_timer('GRO',t1)


!   CALCULATE ROOT-SOIL C AND NUTRIENT EXCHANGE FOR ALL PLANT SPECIES
!   IN 'EXTRACT'
!
  if(lverb)WRITE(*,333)'EXTR'
!    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL EXTRACT(I,J,NHW,NHE,NVN,NVS)
  call end_timer('EXTR',t1)

  end subroutine PlantModel


end module PlantAPI
