module InitNoSaltTransportMod
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: destroy,endrun
  USE MiniMathMod,      ONLY: AZMAX1, fixnegmass, flux_mass_limiter,AZERO
  use TracerPropMod,    only: MolecularWeight
  use EcoSiMParDataMod, only: micpar
  use DebugToolMod
  use SOMDataType
  use ChemTranspDataType
  use GridConsts
  use EcoSIMSolverPar
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use EcosimConst
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use LandSurfDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use PlantDataRateType
  use GridDataType
  use TranspNoSaltDataMod
  use IrrigationDataType
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__
  public :: InitTranspNoSaltModel
  contains
!------------------------------------------------------------------------------------------

  subroutine InitTranspNoSaltModel(I,J,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) :: I, J,NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='InitTranspNoSaltModel'
  integer :: NX,NY

  call PrintInfo('beg '//subname)

  call CopyTracerStates(I,J,NHW,NHE,NVN,NVS)

  call ZeroTracerFluxes(I,J,NHW,NHE,NVN,NVS)

  call SoluteWetDeposition(I,J,NHW,NHE,NVN,NVS)

  call StageSurfaceFluxes(I,J,NY,NX,NHW,NHE,NVN,NVS)

  call StageSoilFluxes(I,J,NHW,NHE,NVN,NVS)

  call PrintInfo('end '//subname)
  end subroutine InitTranspNoSaltModel

!------------------------------------------------------------------------------------------

  subroutine ZeroTracerFluxes(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      trcg_AquaAdv_flxM_snvr(:,1,NY,NX)=0._r8

      trcg_AquaAdv_flx_snvr(idg_beg:idg_NH3,1,NY,NX)          = 0.0_r8
      trcn_AquaAdv_flx_snvr(ids_nut_beg:ids_nuts_end,1,NY,NX) = 0.0_r8

      trcg_Precip2LitrM(idg_beg:idg_NH3,NY,NX)=0._r8
      trcs_Precip2MicpM(ids_beg:ids_end,NY,NX)=0._r8

    ENDDO
  ENDDO      
  end subroutine ZeroTracerFluxes  
!------------------------------------------------------------------------------------------

  subroutine SoluteWetDeposition(I,J,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='SoluteWetDeposition'
  integer :: NY,NX
  integer :: idg

  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

        IF(SnoFalPrec_col(NY,NX).GT.0.0_r8 .OR.   &   !snowfall
          (RainFalPrec_col(NY,NX).GT.0.0_r8 .AND. &   !rainfall and significant snowpack
          VLSnowHeatCapM_snvr(1,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)*1.e-4_r8))THEN !rainfall with snowpack

          call TracerFall2Snowpack(I,J,NY,NX)

          !
          !     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
          !     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
          !     ENTERED IN WEATHER AND IRRIGATION FILES
          !
        ELSEIF((PrecAtm_col(NY,NX).GT.0.0_r8 .OR. IrrigSurface_col(NY,NX).GT.0.0_r8) & !Irrigation or rainfall
          .AND. VLSnowHeatCapM_snvr(1,1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX))THEN     !no snowpack
          !
          !     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
          !     IF SNOWFALL AND IRRIGATION IS ZERO AND SNOWPACK IS ABSENT
          !
          call TracerFall2Grnd(I,J,NY,NX)

          !
          !     NO SOLUTE FLUXES FROM ATMOSPHERE
          !
        ENDIF
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SoluteWetDeposition

!------------------------------------------------------------------------------------------

  subroutine StageSurfaceFluxes(I,J,NY,NX,NHW,NHE,NVN,NVS)
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='StageSurfaceFluxes'
  integer :: NY, NX
  INTEGER :: K,L,ids,idg,idn,idom

  call PrintInfo('beg '//subname)
  !DOM fromp precipiation will be implemented later.
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      DO idg=idg_beg,idg_NH3
        GasDiff2Litr_flx_col(idg,NY,NX)                    = 0._r8
        Gas_WetDepo2Litr_col(idg,NY,NX)                    = 0._r8
        Gas_litr2Soil_flx_col(idg,NY,NX)                   = 0._r8
        GasHydroLoss_litr_flx_col(idg,NY,NX)               = 0._r8
        GasDiff2Soil_flx_col(idg,NY,NX)                    = 0._r8
        Gas_WetDepo2Soil_col(idg,NY,NX)                    = 0._r8
        RGasNetProdSoil_col(idg,NY,NX)                     = 0._r8
        trc_topsoil_flx_col(idg,NY,NX)                     = 0._r8
        TranspNetSoil_flx_col(idg,NY,NX)                   = 0._r8
        TranspNetSoil_flx2_col(idg,NY,NX)                  = 0._r8
        Gas_WetDepo2Snow_col(idg,NY,NX)                    = 0._r8
        Gas_Snowloss_col(idg,NY,NX)                        = 0._r8
        transp_diff_slow_vr(idg,NU(NY,NX):NL(NY,NX),NY,NX) = 0._r8
      ENDDO

      DO  K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end
          DOM_Flo2LitrM(idom,K,NY,NX)    = 0._r8
          DOM_Flo2TopSoilM(idom,K,NY,NX) = 0._r8
        ENDDO
      enddo

      !
      !     GAS AND SOLUTE SINKS AND SOURCES IN SOIL LAYERS FROM MICROBIAL
      !     TRANSFORMATIONS IN 'NITRO' + ROOT EXCHANGE IN 'EXTRACT'
      !     + EQUILIBRIA REACTIONS IN 'SOLUTE' AT SUB-HOURLY TIME STEP
      !
      !     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
      !
      DO ids=ids_beg,ids_end
        SoluteDifusivitytscaledM_vr(ids,0,NY,NX)=SoluteDifusvty_vr(ids,0,NY,NX)*dts_HeatWatTP
      ENDDO

      DO idom=idom_beg,idom_end
        DOM_diffusivitytscaledM_vr(idom,0,NY,NX)=DOMdiffusivity_vr(idom,0,NY,NX)*dts_HeatWatTP
      enddo
      
      DO idg=idg_beg,idg_NH3-1
        if(idg/=idg_O2)RBGCSinkGasMM_vr(idg,0,NY,NX) = trcs_RMicbUptake_vr(idg,0,NY,NX)*dts_gas
      enddo

      RBGCSinkSoluteM_vr(idg_NH3,0,NY,NX)   = -TRChem_sol_NH3_soil_vr(0,NY,NX)*dts_HeatWatTP
      RBGCSinkSoluteM_vr(ids_NH4,0,NY,NX)   = -(RNut_MicbRelease_vr(ids_NH4,0,NY,NX)+trcn_GeoChem_soil_vr(ids_NH4,0,NY,NX))*dts_HeatWatTP
      RBGCSinkSoluteM_vr(ids_NO3,0,NY,NX)   = -(RNut_MicbRelease_vr(ids_NO3,0,NY,NX)+trcn_GeoChem_soil_vr(ids_NO3,0,NY,NX))*dts_HeatWatTP
      RBGCSinkSoluteM_vr(ids_NO2,0,NY,NX)   = -(RNut_MicbRelease_vr(ids_NO2,0,NY,NX)+trcn_GeoChem_soil_vr(ids_NO2,0,NY,NX))*dts_HeatWatTP
      RBGCSinkSoluteM_vr(ids_H2PO4,0,NY,NX) = -(RNut_MicbRelease_vr(ids_H2PO4,0,NY,NX)+trcn_GeoChem_soil_vr(ids_H2PO4,0,NY,NX))*dts_HeatWatTP
      RBGCSinkSoluteM_vr(ids_H1PO4,0,NY,NX) = -(RNut_MicbRelease_vr(ids_H1PO4,0,NY,NX)+trcn_GeoChem_soil_vr(ids_H1PO4,0,NY,NX))*dts_HeatWatTP
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
!
  end subroutine StageSurfaceFluxes


!------------------------------------------------------------------------------------------

  subroutine CopyTracerStates(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='CopyTracerStates'
  integer :: L,K,idg,ids,idom
  integer :: NY,NX

  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      DO L=1,JS
        trcg_solsml2_snvr(idg_beg:idg_NH3,L,NY,NX)          = trcg_solsml_snvr(idg_beg:idg_NH3,L,NY,NX)
        trcn_solsml2_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = trcn_solsml_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX)
      enddo

      ! exclude banded nutrient
      DO idg=idg_beg,idg_end
        trcs_solml2_vr(idg,0,NY,NX)=AZERO(trcs_solml_vr(idg,0,NY,NX))
      ENDDO

    !reset DOM to value before the iteration
      D9979: DO K=1,jcplx
        DO idom=idom_beg,idom_end
          DOM_MicP2_vr(idom,K,0,NY,NX)=AZERO(DOM_MicP_vr(idom,K,0,NY,NX))
        ENDDO
      ENDDO D9979

      DO L=NU(NY,NX),NL(NY,NX)

        !
        !     STATE VARIABLES FOR GASES AND SOLUTES USED IN 'TranspNoSalt'
        !     TO STORE SUB-HOURLY CHANGES DURING FLUX CALCULATIONS
        !     INCLUDING TRANSFORMATIONS FROM NITRO.F, UPTAKE.F AND SOLUTE.F
        ! exclude NH3B
        DO idg=idg_beg,idg_NH3
          trcg_gasml2_vr(idg,L,NY,NX)=AZERO(trcg_gasml_vr(idg,L,NY,NX))
        ENDDO

        !reset DOM to value before the iteration
        !Concern: RDOMMicProd_vr includes microbial modification of DOM
        !it does not include exchange among different complexes (as one kind of microbial priming). 

        DO  K=1,jcplx
          DO idom=idom_beg,idom_end
            DOM_MicP2_vr(idom,K,L,NY,NX) = AZERO(DOM_MicP_vr(idom,K,L,NY,NX))
            DOM_MacP2_vr(idom,K,L,NY,NX)    = AZERO(DOM_MacP_vr(idom,K,L,NY,NX))
          ENDDO
        enddo

        DO ids=ids_beg,ids_end
          trcs_solml2_vr(ids,L,NY,NX) = AZERO(trcs_solml_vr(ids,L,NY,NX))    !micropore solute tracer
          trcs_soHml2_vr(ids,L,NY,NX) = AZERO(trcs_soHml_vr(ids,L,NY,NX))    !macropore solute tracer

          if(trcs_solml2_vr(ids,L,NY,NX)<0._r8)then
            write(*,*)I*1000+J,'micv2 ',trcs_names(ids),L,trcs_solml2_vr(ids,L,NY,NX),trcs_solml_vr(ids,L,NY,NX)
            call endrun(trim(mod_filename)//' at line',__LINE__)        
          endif
          if(trcs_soHml_vr(ids,L,NY,NX)<0._r8)then
            write(*,*)'macv2 ',trcs_names(ids),L,trcs_soHml_vr(ids,L,NY,NX)
            call endrun(trim(mod_filename)//' at line',__LINE__)                
          endif
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine CopyTracerStates
!------------------------------------------------------------------------------------------
  subroutine StageSoilFluxes(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='StageSoilFluxes'
  integer :: NY,NX
  integer :: L,K,idg,ids,idom

  !dts_gas=1/NPG with NPG=NPH*NPT
  !dts_HeatWatTP=1/NPT  

  call PrintInfo('beg '//subname)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU(NY,NX),NL(NY,NX)        

        DO idg=idg_beg,idg_NH3-1
          if(idg/=idg_O2)then      
            RBGCSinkGasMM_vr(idg,L,NY,NX) = (trcs_RMicbUptake_vr(idg,L,NY,NX)-trcs_deadroot2soil_vr(idg,L,NY,NX))*dts_gas
          endif
        enddo

        RBGCSinkGasMM_vr(idg_CO2,L,NY,NX) = RBGCSinkGasMM_vr(idg_CO2,L,NY,NX)- &
          (TProd_CO2_geochem_soil_vr(L,NY,NX)+RootCO2Ar2Soil_vr(L,NY,NX))*dts_gas

        RBGCSinkGasMM_vr(idg_N2,L,NY,NX)     = RBGCSinkGasMM_vr(idg_N2,L,NY,NX)+RootN2Fix_vr(L,NY,NX)*dts_gas
        RBGCSinkGasMM_vr(idg_NH3,L,NY,NX)    = -TRChem_gas_NH3_geochem_vr(L,NY,NX)*dts_gas                    !geochemical NH3 source

        RBGCSinkSoluteM_vr(idg_NH3B,L,NY,NX) = (-trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX))*dts_HeatWatTP
        RBGCSinkSoluteM_vr(idg_NH3,L,NY,NX)  = (-TRChem_sol_NH3_soil_vr(L,NY,NX)-trcs_deadroot2soil_vr(idg_NH3,L,NY,NX))*dts_HeatWatTP
        DO ids=ids_nut_beg,ids_nuts_end
          RBGCSinkSoluteM_vr(ids,L,NY,NX) = (-RNut_MicbRelease_vr(ids,L,NY,NX)-trcn_GeoChem_soil_vr(ids,L,NY,NX)+trcs_plant_uptake_vr(ids,L,NY,NX))*dts_HeatWatTP
        ENDDO

        DO ids=ids_NH4B,ids_nutb_end
          RBGCSinkSoluteM_vr(ids,L,NY,NX) = (-RNut_MicbRelease_vr(ids,L,NY,NX)-trcn_RChem_band_soil_vr(ids,L,NY,NX)+trcs_plant_uptake_vr(ids,L,NY,NX))*dts_HeatWatTP
        ENDDO
        !
        !     SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
        !
        do idg=idg_beg,idg_NH3
          trcs_Irrig_flx_vr(idg,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*trcg_irrig_mole_conc_col(idg,NY,NX)*MolecularWeight(idg)
        ENDDO
    
        trcs_Irrig_flx_vr(ids_NH4,L,NY,NX)    = FWatIrrigate2MicP_vr(L,NY,NX)*NH4_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
        trcs_Irrig_flx_vr(idg_NH3,L,NY,NX)    = FWatIrrigate2MicP_vr(L,NY,NX)*trcg_irrig_mole_conc_col(idg_NH3,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
        trcs_Irrig_flx_vr(ids_NO3,L,NY,NX)    = FWatIrrigate2MicP_vr(L,NY,NX)*NO3_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)*natomw
        trcs_Irrig_flx_vr(ids_H1PO4,L,NY,NX)  = FWatIrrigate2MicP_vr(L,NY,NX)*HPO4_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
        trcs_Irrig_flx_vr(ids_H2PO4,L,NY,NX)  = FWatIrrigate2MicP_vr(L,NY,NX)*H2PO4_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
        trcs_Irrig_flx_vr(ids_NH4B,L,NY,NX)   = FWatIrrigate2MicP_vr(L,NY,NX)*NH4_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
        trcs_Irrig_flx_vr(idg_NH3B,L,NY,NX)   = FWatIrrigate2MicP_vr(L,NY,NX)*trcg_irrig_mole_conc_col(idg_NH3,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
        trcs_Irrig_flx_vr(ids_NO3B,L,NY,NX)   = FWatIrrigate2MicP_vr(L,NY,NX)*NO3_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)*natomw
        trcs_Irrig_flx_vr(ids_H1PO4B,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*HPO4_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw
        trcs_Irrig_flx_vr(ids_H2PO4B,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*H2PO4_irrig_mole_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw
        !
        !     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
        !
        DO ids=ids_beg,ids_end
          trcsol_Irrig_flxM_vr(ids,L,NY,NX)=trcs_Irrig_flx_vr(ids,L,NY,NX)*dts_HeatWatTP
        ENDDO
        !
        !     GAS AND SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
        !
        DO idom=idom_beg,idom_end
          DOM_diffusivitytscaledM_vr(idom,L,NY,NX)=DOMdiffusivity_vr(idom,L,NY,NX)*dts_HeatWatTP
        enddo

        DO ids=ids_beg,ids_end
          SoluteDifusivitytscaledM_vr(ids,L,NY,NX)=SoluteDifusvty_vr(ids,L,NY,NX)*dts_HeatWatTP
        ENDDO

        DO idg=idg_beg,idg_NH3
          GasDifctScaledMM_vr(idg,L,NY,NX)=GasDifc_vr(idg,L,NY,NX)*dts_gas
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine StageSoilFluxes

!------------------------------------------------------------------------------------------

  subroutine TracerFall2Snowpack(I,J,NY,NX)
  implicit none

  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  integer :: idg

  DO idg=idg_beg,idg_NH3
    trcg_AquaAdv_flxM_snvr(idg,1,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcg_rain_mole_conc_col(idg,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcg_irrig_mole_conc_col(idg,NY,NX))*MolecularWeight(idg)*dts_HeatWatTP
  ENDDO

  trcn_AquaAdv_flxM_snvr(ids_NH4,1,NY,NX)   = (Rain2SoilSurf_col(NY,NX)*NH4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NH4_irrig_mole_conc(I,NY,NX))*natomw*dts_HeatWatTP
  trcn_AquaAdv_flxM_snvr(ids_NO3,1,NY,NX)   = (Rain2SoilSurf_col(NY,NX)*NO3_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NO3_irrig_mole_conc(I,NY,NX))*natomw*dts_HeatWatTP
  trcn_AquaAdv_flxM_snvr(ids_H1PO4,1,NY,NX) = (Rain2SoilSurf_col(NY,NX)*HPO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*HPO4_irrig_mole_conc(I,NY,NX))*patomw*dts_HeatWatTP
  trcn_AquaAdv_flxM_snvr(ids_H2PO4,1,NY,NX) = (Rain2SoilSurf_col(NY,NX)*H2PO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*H2PO4_irrig_mole_conc(I,NY,NX))*patomw*dts_HeatWatTP

  !     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
  !     IF RAINFALL AND IRRIGATION IS ZERO IF SNOWPACK IS PRESENT
  !

  end subroutine TracerFall2Snowpack
!------------------------------------------------------------------------------------------

  subroutine TracerFall2Grnd(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer :: idg

  trcg_AquaAdv_flxM_snvr(idg_beg:idg_NH3,1,NY,NX)          = 0.0_r8
  trcn_AquaAdv_flxM_snvr(ids_nut_beg:ids_nuts_end,1,NY,NX) = 0.0_r8
    
  do idg=idg_beg,idg_NH3
    trcg_Precip2LitrM(idg,NY,NX)=dts_HeatWatTP*(Rain2LitRSurf_col(NY,NX)*trcg_rain_mole_conc_col(idg,NY,NX) &
      +Irrig2LitRSurf_col(NY,NX)*trcg_irrig_mole_conc_col(idg,NY,NX))*MolecularWeight(idg)

    trcs_Precip2MicpM(idg,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*trcg_rain_mole_conc_col(idg,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcg_irrig_mole_conc_col(idg,NY,NX))*MolecularWeight(idg)
  ENDDO

  trcg_Precip2LitrM(ids_NH4,NY,NX)=dts_HeatWatTP*(Rain2LitRSurf_col(NY,NX)*NH4_rain_mole_conc(NY,NX) &
    +Irrig2LitRSurf_col(NY,NX)*NH4_irrig_mole_conc(I,NY,NX))*natomw

  trcg_Precip2LitrM(idg_NH3,NY,NX)=dts_HeatWatTP*(Rain2LitRSurf_col(NY,NX)*trcg_rain_mole_conc_col(idg_NH3,NY,NX) &
    +Irrig2LitRSurf_col(NY,NX)*trcg_irrig_mole_conc_col(idg_NH3,NY,NX))*natomw

  trcg_Precip2LitrM(ids_NO3,NY,NX)=dts_HeatWatTP*(Rain2LitRSurf_col(NY,NX)*NO3_rain_mole_conc(NY,NX) &
    +Irrig2LitRSurf_col(NY,NX)*NO3_irrig_mole_conc(I,NY,NX))*natomw

  trcg_Precip2LitrM(ids_NO2,NY,NX)=0.0_r8

  trcg_Precip2LitrM(ids_H1PO4,NY,NX)=dts_HeatWatTP*(Rain2LitRSurf_col(NY,NX)*HPO4_rain_mole_conc(NY,NX) &
    +Irrig2LitRSurf_col(NY,NX)*HPO4_irrig_mole_conc(I,NY,NX))*patomw

  trcg_Precip2LitrM(ids_H2PO4,NY,NX)=dts_HeatWatTP*(Rain2LitRSurf_col(NY,NX)*H2PO4_rain_mole_conc(NY,NX) &
    +Irrig2LitRSurf_col(NY,NX)*H2PO4_irrig_mole_conc(I,NY,NX))*patomw

  trcs_Precip2MicpM(ids_NH4,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*NH4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NH4_irrig_mole_conc(I,NY,NX))*natomw*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(idg_NH3,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*trcg_rain_mole_conc_col(idg_NH3,NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*trcg_irrig_mole_conc_col(idg_NH3,NY,NX))*natomw*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(ids_NO3,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*NO3_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NO3_irrig_mole_conc(I,NY,NX))*natomw*trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(ids_NO2,NY,NX)=0.0_r8

  trcs_Precip2MicpM(ids_H1PO4,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*HPO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*HPO4_irrig_mole_conc(I,NY,NX))*patomw*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(ids_H2PO4,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*H2PO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*H2PO4_irrig_mole_conc(I,NY,NX))*patomw*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(ids_NH4B,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*NH4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NH4_irrig_mole_conc(I,NY,NX))*natomw*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(idg_NH3B,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*trcg_rain_mole_conc_col(idg_NH3,NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*trcg_irrig_mole_conc_col(idg_NH3,NY,NX))*natomw*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(ids_NO3B,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*NO3_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NO3_irrig_mole_conc(I,NY,NX))*natomw*trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(ids_NO2B,NY,NX)=0.0_r8

  trcs_Precip2MicpM(ids_H1PO4B,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*HPO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*HPO4_irrig_mole_conc(I,NY,NX))*patomw*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)

  trcs_Precip2MicpM(ids_H2PO4B,NY,NX)=dts_HeatWatTP*(Rain2SoilSurf_col(NY,NX)*H2PO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*H2PO4_irrig_mole_conc(I,NY,NX))*patomw*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)

  end subroutine TracerFall2Grnd            

end module InitNoSaltTransportMod  