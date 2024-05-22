  module YearMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use EcoSIMCtrlMod
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use PlantTraitDataType
  use SurfLitterDataType, only : XCORP
  use PlantDataRateType
  use CanopyDataType
  use RootDataType
  use PlantMgmtDataType
  use SOMDataType
  use EcosimBGCFluxType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use SedimentDataType
  use GridDataType
  use EcoSIMConfig  
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: SetAnnualAccumlators
  contains
!------------------------------------------------------------------------------------------

  subroutine SetAnnualAccumlators(NY,NX)
  
  !set accumulators at the beginning of year
  implicit none
  integer, intent(in) :: NY,NX

  integer :: NZ,NE

  AmendCFlx_col(NY,NX)=0._r8
  LiterfalOrgM_col(ielmc,NY,NX)=0._r8
  UCOP(NY,NX)=0._r8
  Eco_NBP_col(NY,NX)=0._r8
  URAIN_col(NY,NX)=0._r8
  UEVAP(NY,NX)=0._r8
  URUN(NY,NX)=0._r8
  USEDOU(NY,NX)=0._r8
  AnualH2OLoss_col(NY,NX)=0._r8
  HydroIonFlx_col(NY,NX)=0._r8
  FertNFlx_col(NY,NX)=0._r8
  LiterfalOrgM_col(ielmn,NY,NX)=0._r8
  HydroSufDONFlx_col(NY,NX)=0._r8
  HydroSufDINFlx_col(NY,NX)=0._r8
  FerPFlx_col(NY,NX)=0._r8
  LiterfalOrgM_col(ielmp,NY,NX)=0._r8
  HydroSufDOPFlx_col(NY,NX)=0._r8
  HydroSufDIPFlx_col(NY,NX)=0._r8
  CO2byFire_col(NY,NX)=0._r8
  CH4byFire_col(NY,NX)=0._r8
  O2byFire_col(NY,NX)=0._r8
  N2ObyFire_col(NY,NX)=0._r8
  NH3byFire_col(NY,NX)=0._r8
  PO4byFire_col(NY,NX)=0._r8
  Eco_HR_col(NY,NX)=0._r8
  Eco_GPP_col(NY,NX)=0._r8
  Eco_NPP_col(NY,NX)=0._r8
  Eco_AutoR_col(NY,NX)=0._r8
  EcoHavstElmnt_col(:,NY,NX)=0._r8
  NetNH4Mineralize_col(NY,NX)=0._r8
  NetPO4Mineralize_col(NY,NX)=0._r8

  D960: DO NZ=1,NP0(NY,NX)
  !NetCumElmntFlx2Plant_pft: effect of canopy element status on seed set
    DO NE=1,NumPlantChemElms
      NetCumElmntFlx2Plant_pft(NE,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(NE,NZ,NY,NX) &
        +PlantExudChemElmCum_pft(NE,NZ,NY,NX)-LitrfalStrutElmsCum_pft(NE,NZ,NY,NX)
      EcoHavstElmntCum_pft(NE,NZ,NY,NX)=EcoHavstElmntCum_pft(NE,NZ,NY,NX)+EcoHavstElmnt_pft(NE,NZ,NY,NX)
      NodulInfectElmsCum_pft(NE,NZ,NY,NX)=0._r8
    ENDDO
    NetCumElmntFlx2Plant_pft(ielmc,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmc,NZ,NY,NX) &
      +GrossCO2FixCum_pft(NZ,NY,NX)+GrossRespCum_pft(NZ,NY,NX) &
      -CO2ByFire_pft(NZ,NY,NX)-CH4ByFire_pft(NZ,NY,NX)
    NetCumElmntFlx2Plant_pft(ielmn,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmn,NZ,NY,NX) &
      +NH3EmisCum_pft(NZ,NY,NX)+PlantN2FixCum_pft(NZ,NY,NX) &
      -NH3byFire_pft(NZ,NY,NX)-N2ObyFire_pft(NZ,NY,NX)
    NetCumElmntFlx2Plant_pft(ielmp,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmp,NZ,NY,NX) &
      -PO4byFire_pft(NZ,NY,NX)

! the following variables are accumulated daily
    GrossCO2FixCum_pft(NZ,NY,NX)=0._r8
    PlantExudChemElmCum_pft(:,NZ,NY,NX)=0._r8
    GrossRespCum_pft(NZ,NY,NX)=0._r8
    CanopyRespC_pft(NZ,NY,NX)=0._r8
    ETCanopy_pft(NZ,NY,NX)=0._r8
    PlantN2FixCum_pft(NZ,NY,NX)=0._r8
    NH3EmisCum_pft(NZ,NY,NX)=0._r8
    CO2ByFire_pft(NZ,NY,NX)=0._r8
    CH4ByFire_pft(NZ,NY,NX)=0._r8
    O2ByFire_pft(NZ,NY,NX)=0._r8
    NH3byFire_pft(NZ,NY,NX)=0._r8
    N2ObyFire_pft(NZ,NY,NX)=0._r8
    PO4byFire_pft(NZ,NY,NX)=0._r8
    EcoHavstElmnt_pft(:,NZ,NY,NX)=0._r8
    SurfLitrfalStrutElmsCum_pft(:,NZ,NY,NX)=0._r8
    LitrfalStrutElmsCum_pft(:,NZ,NY,NX)=0._r8
  ENDDO D960
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    TSED(NY,NX)=0._r8
  ENDIF
  end subroutine SetAnnualAccumlators        

end module YearMod