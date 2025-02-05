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

  subroutine SetAnnualAccumlators(I, NHW, NHE, NVN, NVS)
  
  !set accumulators at the beginning of year
  implicit none
  integer, intent(in) :: I, NHW, NHE, NVN, NVS
  integer :: NY,NX

  integer :: NZ,NE

  D955: DO NX=NHW,NHE
    D950: DO NY=NVN,NVS
!     RESET ANNUAL FLUX ACCUMULATORS AT START OF ANNUAL CYCLE
!     ALAT=latitude +ve=N,-ve=S
!

      IF((ALAT(NY,NX).GE.0.0_r8.AND.I.EQ.1) .OR. (ALAT(NY,NX).LT.0.0_r8.AND.I.EQ.1))THEN        
        !to be revised for GDD 
        Gas_Prod_TP_cumRes_col(:,NY,NX)     = 0._r8
        QdewCanopy_CumYr_pft(:,NY,NX)     = 0._r8
        GDD_col(NY,NX)                    = 0._r8
        AmendCFlx_CumYr_col(NY,NX)        = 0._r8
        LiterfalOrgM_col(:,NY,NX)         = 0._r8
        RootResp_CumYr_col(NY,NX)         = 0._r8
        Eco_NBP_CumYr_col(NY,NX)          = 0._r8
        QRain_CumYr_col(NY,NX)            = 0._r8
        QEvap_CumYr_col(NY,NX)            = 0._r8
        Qrunoff_CumYr_col(NY,NX)          = 0._r8
        SedmErossLoss_CumYr_col(NY,NX)    = 0._r8
        H2OLoss_CumYr_col(NY,NX)          = 0._r8
        HydroIonFlx_CumYr_col(NY,NX)      = 0._r8
        FertNFlx_CumYr_col(NY,NX)         = 0._r8
        HydroSufDONFlx_CumYr_col(NY,NX)   = 0._r8
        HydroSufDINFlx_CumYr_col(NY,NX)   = 0._r8
        FerPFlx_CumYr_col(NY,NX)          = 0._r8
        HydroSufDOPFlx_CumYr_col(NY,NX)   = 0._r8
        HydroSufDIPFlx_CumYr_col(NY,NX)   = 0._r8
        CO2byFire_CumYr_col(NY,NX)        = 0._r8
        CH4byFire_CumYr_col(NY,NX)        = 0._r8
        O2byFire_CumYr_col(NY,NX)         = 0._r8
        N2ObyFire_CumYr_col(NY,NX)        = 0._r8
        NH3byFire_CumYr_col(NY,NX)        = 0._r8
        PO4byFire_CumYr_col(NY,NX)        = 0._r8
        Eco_HR_CumYr_col(NY,NX)           = 0._r8
        Eco_GPP_CumYr_col(NY,NX)          = 0._r8
        Eco_NPP_CumYr_col(NY,NX)          = 0._r8
        Eco_AutoR_CumYr_col(NY,NX)        = 0._r8
        EcoHavstElmnt_CumYr_col(:,NY,NX)  = 0._r8
        NetNH4Mineralize_CumYr_col(NY,NX) = 0._r8
        NetPO4Mineralize_CumYr_col(NY,NX) = 0._r8
        HoursTooLowPsiCan_pft(:,NY,NX)    = 0._r8
        D960: DO NZ=1,NP0(NY,NX)
          !NetCumElmntFlx2Plant_pft: effect of canopy element status on seed set
          DO NE=1,NumPlantChemElms
            NetCumElmntFlx2Plant_pft(NE,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(NE,NZ,NY,NX) &
              +PlantExudElm_CumYr_pft(NE,NZ,NY,NX)-LitrfalStrutElms_CumYr_pft(NE,NZ,NY,NX)
              
            EcoHavstElmntCum_pft(NE,NZ,NY,NX)=EcoHavstElmntCum_pft(NE,NZ,NY,NX)+EcoHavstElmnt_CumYr_pft(NE,NZ,NY,NX)
            NodulInfectElmsCum_pft(NE,NZ,NY,NX)=0._r8
          ENDDO
          NetCumElmntFlx2Plant_pft(ielmc,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmc,NZ,NY,NX) &
            +GrossCO2Fix_CumYr_pft(NZ,NY,NX)+GrossRespC_CumYr_pft(NZ,NY,NX) &
            -CO2ByFire_CumYr_pft(NZ,NY,NX)-CH4ByFire_CumYr_pft(NZ,NY,NX)
          NetCumElmntFlx2Plant_pft(ielmn,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmn,NZ,NY,NX) &
            +NH3Emis_CumYr_pft(NZ,NY,NX)+PlantN2Fix_CumYr_pft(NZ,NY,NX) &
            -NH3byFire_CumYr_pft(NZ,NY,NX)-N2ObyFire_CumYr_pft(NZ,NY,NX) &
            +RootUptk_N_CumYr_pft(NZ,NY,NX)
          NetCumElmntFlx2Plant_pft(ielmp,NZ,NY,NX)=NetCumElmntFlx2Plant_pft(ielmp,NZ,NY,NX) &
            -PO4byFire_CumYr_pft(NZ,NY,NX)+RootUptk_P_CumYr_pft(NZ,NY,NX)

          RootUptk_N_CumYr_pft(NZ,NY,NX)             = 0._r8
          RootUptk_P_CumYr_pft(NZ,NY,NX)             = 0._r8
          GrossCO2Fix_CumYr_pft(NZ,NY,NX)            = 0._r8
          PlantExudElm_CumYr_pft(:,NZ,NY,NX)         = 0._r8
          GrossRespC_CumYr_pft(NZ,NY,NX)             = 0._r8
          CanopyRespC_CumYr_pft(NZ,NY,NX)            = 0._r8
          ETCanopy_CumYr_pft(NZ,NY,NX)               = 0._r8
          PlantN2Fix_CumYr_pft(NZ,NY,NX)             = 0._r8
          NH3Emis_CumYr_pft(NZ,NY,NX)                = 0._r8
          CO2ByFire_CumYr_pft(NZ,NY,NX)              = 0._r8
          CH4ByFire_CumYr_pft(NZ,NY,NX)              = 0._r8
          O2ByFire_CumYr_pft(NZ,NY,NX)               = 0._r8
          NH3byFire_CumYr_pft(NZ,NY,NX)              = 0._r8
          N2ObyFire_CumYr_pft(NZ,NY,NX)              = 0._r8
          PO4byFire_CumYr_pft(NZ,NY,NX)              = 0._r8
          EcoHavstElmnt_CumYr_pft(:,NZ,NY,NX)        = 0._r8
          SurfLitrfalStrutElms_CumYr_pft(:,NZ,NY,NX) = 0._r8
          LitrfalStrutElms_CumYr_pft(:,NZ,NY,NX)     = 0._r8
        ENDDO D960
        IF(iErosionMode.EQ.ieros_frzthaweros .OR. iErosionMode.EQ.ieros_frzthawsomeros)THEN
          TSED(NY,NX)=0._r8
        ENDIF
      ENDIF
    ENDDO D950
  ENDDO D955

  end subroutine SetAnnualAccumlators        

end module YearMod