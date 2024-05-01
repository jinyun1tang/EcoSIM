module MicBGCAPI

  use data_kind_mod  , only : r8 => DAT_KIND_R8
  use NitrosMod      , only : VerticalLitterMixLvsLL
  use NitroDisturbMod, only : SOMRemovalByDisturbance
  use EcoSIMSolverPar
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod , only : micforctype
  use minimathmod    , only : AZMAX1
  use TracerIDMod
  use SoilBGCDataType
  USE PlantDataRateType
  USE SoilWaterDataType
  USE SurfLitterDataType
  USE EcosimBGCFluxType
  USE EcoSIMCtrlDataType
  USE ClimForcDataType
  USE GridDataType
  USE SoilPropertyDataType
  use SoilPhysDataType
  use SOMDataType
  use ChemTranspDataType
  use MicrobialDataType
  use IrrigationDataType
  use SoilHeatDataType
  use EcoSiMParDataMod, only : micpar
  use MicBGCMod, only : SoilBGCOneLayer
implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  type(micforctype) :: micfor
  type(micsttype) :: micstt
  type(micfluxtype) :: micflx
  integer :: curI,curJ


  public :: MicrobeModel
  public :: MicAPI_Init
  public :: MicAPI_cleanup
  contains

!------------------------------------------------------------------------------------------

  subroutine MicAPI_Init()

  implicit none

  call micfor%Init()
  call micstt%Init()
  call micflx%Init()

  end subroutine MicAPI_Init
!------------------------------------------------------------------------------------------
  subroutine MicAPI_cleanup()

  implicit none


  call micfor%destroy()
  call micflx%destroy()
  call micstt%destroy()
  end subroutine MicAPI_cleanup
!------------------------------------------------------------------------------------------
  subroutine MicrobeModel(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOIL BIOLOGICAL TRANSFORMATIONS
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: L,NX,NY

!   begin_execution
  curI=I; curJ=J
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
!
!       VOLWZ=water volume used to calculate aqueous microbial
!       concentrations that drive microbial density effects on
!       decomposition
      D998: DO L=0,NL(NY,NX)
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN
             call MicBGC1Layer(I,J,L,NY,NX)
          ELSE
            trcg_RMicbTransf_vr(idg_CO2,L,NY,NX)=0.0_r8
            trcg_RMicbTransf_vr(idg_CH4,L,NY,NX)=0.0_r8
            trcg_RMicbTransf_vr(idg_H2,L,NY,NX)=0.0_r8
            trcg_RMicbTransf_vr(idg_O2,L,NY,NX)=0.0_r8
            trcg_RMicbTransf_vr(idg_N2,L,NY,NX)=0.0_r8
            trcg_RMicbTransf_vr(idg_N2O,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_NH4,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_NO3,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_NO2,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_H2PO4,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_H1PO4,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_NH4B,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_NO3B,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_NO2B,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX)=0.0_r8
            RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX)=0.0_r8
            Micb_N2Fixation_vr(L,NY,NX)=0.0_r8
          ENDIF
  !     MIX LITTER C BETWEEN ADJACENT SOIL LAYERS L AND LL
          call VerticalLitterMixLvsLL(I,J,L,NY,NX)

        ELSE
          trcg_RMicbTransf_vr(idg_CO2,L,NY,NX)=0.0_r8
          trcg_RMicbTransf_vr(idg_CH4,L,NY,NX)=0.0_r8
          trcg_RMicbTransf_vr(idg_H2,L,NY,NX)=0.0_r8
          trcg_RMicbTransf_vr(idg_O2,L,NY,NX)=0.0_r8
          trcg_RMicbTransf_vr(idg_N2,L,NY,NX)=0.0_r8
          trcg_RMicbTransf_vr(idg_N2O,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_NH4,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_NO3,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_NO2,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_H2PO4,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_H1PO4,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_NH4B,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_NO3B,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_NO2B,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX)=0.0_r8
          RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX)=0.0_r8
          Micb_N2Fixation_vr(L,NY,NX)=0.0_r8
        ENDIF
      ENDDO D998
!
!       SOC LOSS IF FIRE OR REMOVAL EVENT IS ENTERED IN DISTURBANCE FILE
!
      call SOMRemovalByDisturbance(I,J,NY,NX)
    ENDDO D9990
  ENDDO D9995
  RETURN
  END subroutine MicrobeModel

!------------------------------------------------------------------------------------------

  subroutine MicBGC1Layer(I,J,L,NY,NX)

  implicit none
  integer, intent(in) :: I,J,L,NY,NX

  call micflx%ZeroOut()

  call MicAPISend(L,NY,NX,micfor,micstt,micflx)

  call SoilBGCOneLayer(micfor,micstt,micflx)

  call MicAPIRecv(L,NY,NX,micfor%litrm,micstt,micflx)

  end subroutine Micbgc1Layer
!------------------------------------------------------------------------------------------

  subroutine MicAPISend(L,NY,NX,micfor,micstt,micflx)
  implicit none
  integer, intent(in) :: L,NY,NX
  type(micforctype), intent(inout) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx

  integer :: NumMicbFunGroups, jcplx, k_POM, k_humus
  integer :: kk, ndbiomcp, nlbiomcp, NumMicrobAutrophCmplx, NumHetetrMicCmplx
  NumMicbFunGroups=micpar%NumMicbFunGroups
  jcplx=micpar%jcplx

  ndbiomcp = micpar%ndbiomcp
  nlbiomcp = micpar%nlbiomcp
  NumMicrobAutrophCmplx  = micpar%NumMicrobAutrophCmplx
  NumHetetrMicCmplx  = micpar%NumHetetrMicCmplx
  k_humus  = micpar%k_humus
  k_POM    = micpar%k_POM

  micfor%ZERO  =ZERO
  micfor%CCH4E =AtmGasCgperm3(idg_CH4,NY,NX)
  micfor%COXYE =AtmGasCgperm3(idg_O2,NY,NX)
  micfor%O2_irrig_conc  =O2_irrig_conc(NY,NX)
  micfor%O2_rain_conc  =O2_rain_conc(NY,NX)
  micfor%Irrig2LitRSurf =Irrig2LitRSurf(NY,NX)
  micfor%Rain2LitRSurf_col =Rain2LitRSurf_col(NY,NX)
  micfor%TempOffset=TempOffset_col(NY,NX)
  micfor%VLitR  =VLitR(NY,NX)
  micfor%VWatLitRHoldCapcity=VWatLitRHoldCapcity(NY,NX)
  micfor%VWatLitRHoldCapcity=VWatLitRHoldCapcity(NY,NX)
  micfor%ZEROS2=ZEROS2(NY,NX)
  micfor%ZEROS =ZEROS(NY,NX)
  micfor%VLSoilMicP  =VLSoilMicP(L,NY,NX)
  micfor%THETY =THETY(L,NY,NX)
  micfor%POROS =POROS(L,NY,NX)
  micfor%FieldCapacity    =FieldCapacity(L,NY,NX)
  micfor%TKS   =TKS(L,NY,NX)
  micfor%THETW =THETW(L,NY,NX)
  micfor%PH    =PH(L,NY,NX)
  micfor%SoilMicPMassLayer  =SoilMicPMassLayer(L,NY,NX)
  micfor%VLSoilPoreMicP  =VLSoilPoreMicP_vr(L,NY,NX)
  micfor%TScal4Difsvity  =TScal4Difsvity_vr(L,NY,NX)
  micfor%VLNOB =trcs_VLN_vr(ids_NO3B,L,NY,NX)
  micfor%VLNO3 =trcs_VLN_vr(ids_NO3,L,NY,NX)
  micfor%VLNH4 =trcs_VLN_vr(ids_NH4,L,NY,NX)
  micfor%VLNHB =trcs_VLN_vr(ids_NH4B,L,NY,NX)
  micfor%VLPO4 =trcs_VLN_vr(ids_H1PO4,L,NY,NX)
  micfor%VLPOB =trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
  micfor%PSISoilMatricP =PSISoilMatricP(L,NY,NX)
  micfor%O2AquaDiffusvity =SolDifc_vr(idg_O2,L,NY,NX)
  micfor%ORGC  =ORGC(L,NY,NX)
  micfor%RNO2EcoUptkSoilPrev =RNO2EcoUptkSoilPrev_vr(L,NY,NX)
  micfor%RN2OEcoUptkSoilPrev =RN2OEcoUptkSoilPrev_vr(L,NY,NX)
  micfor%RNO2EcoUptkBandPrev =RNO2EcoUptkBandPrev_vr(L,NY,NX)
  micfor%RO2EcoDmndPrev =RO2EcoDmndPrev_vr(L,NY,NX)
  micfor%RO2GasXchangePrev =RO2GasXchangePrev_vr(L,NY,NX)
  micfor%RCH4L =RCH4L(L,NY,NX)
  micfor%RNH4EcoDmndBandPrev =RNH4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RNO3EcoDmndBandPrev =RNO3EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RH2PO4EcoDmndBandPrev =RH2PO4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RH1PO4EcoDmndBandPrev =RH1PO4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RO2AquaXchangePrev =RO2AquaXchangePrev_vr(L,NY,NX)
  micfor%RDOMEcoDmndPrev(1:jcplx)=RDOMEcoDmndPrev_vr(1:jcplx,L,NY,NX)
  micfor%RAcetateEcoDmndPrev(1:jcplx)=RAcetateEcoDmndPrev_vr(1:jcplx,L,NY,NX)
  micfor%litrm=(L==0)
  micfor%Lsurf=(L==NU(NY,NX))
  if(micfor%litrm)then
    micstt%ZNH4TU=AZMAX1(trc_solml_vr(ids_NH4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_NH4B,NU(NY,NX),NY,NX))
    micstt%ZNO3TU=AZMAX1(trc_solml_vr(ids_NO3,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_NO3B,NU(NY,NX),NY,NX))
    micstt%H1P4TU=AZMAX1(trc_solml_vr(ids_H1PO4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_H1PO4B,NU(NY,NX),NY,NX))
    micstt%H2P4TU=AZMAX1(trc_solml_vr(ids_H2PO4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_H2PO4B,NU(NY,NX),NY,NX))
    micstt%CNH4BU=trc_solcl_vr(ids_NH4B,NU(NY,NX),NY,NX)
    micstt%CNH4SU=trc_solcl_vr(ids_NH4,NU(NY,NX),NY,NX)
    micstt%CH2P4U=trc_solcl_vr(ids_H2PO4,NU(NY,NX),NY,NX)
    micstt%CH2P4BU=trc_solcl_vr(ids_H2PO4B,NU(NY,NX),NY,NX)
    micstt%CH1P4U=trc_solcl_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    micstt%CH1P4BU=trc_solcl_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
    micstt%CNO3SU=trc_solcl_vr(ids_NO3,NU(NY,NX),NY,NX)
    micstt%CNO3BU=trc_solcl_vr(ids_NO3B,NU(NY,NX),NY,NX)
    micstt%OSC13U=SolidOM_vr(ielmc,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)
    micstt%OSN13U=SolidOM_vr(ielmn,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)
    micstt%OSP13U=SolidOM_vr(ielmp,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)
    micstt%OSC14U=SolidOM_vr(ielmc,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSN14U=SolidOM_vr(ielmn,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSP14U=SolidOM_vr(ielmp,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSC24U=SolidOM_vr(ielmc,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSN24U=SolidOM_vr(ielmn,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)
    micstt%OSP24U=SolidOM_vr(ielmp,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)
    micfor%RNH4EcoDmndLitrPrev =RNH4EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%RNO3EcoDmndLitrPrev =RNO3EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%RH2PO4EcoDmndLitrPrev =RH2PO4EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%RH1PO4EcoDmndLitrPrev =RH1PO4EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%VOLWU =VLWatMicP(NU(NY,NX),NY,NX)
    micfor%CFOMCU=CFOMC(1:2,NU(NY,NX),NY,NX)
  else
    micfor%CFOMC =CFOMC(1:2,L,NY,NX)
  endif
  micstt%CNH4B =trc_solcl_vr(ids_NH4B,L,NY,NX)
  micstt%CNH4S =trc_solcl_vr(ids_NH4,L,NY,NX)
  micstt%CH2P4 =trc_solcl_vr(ids_H2PO4,L,NY,NX)
  micstt%CH2P4B=trc_solcl_vr(ids_H2PO4B,L,NY,NX)
  micstt%CH1P4=trc_solcl_vr(ids_H1PO4,L,NY,NX)
  micstt%CH1P4B=trc_solcl_vr(ids_H1PO4B,L,NY,NX)
  micstt%CNO3S=trc_solcl_vr(ids_NO3,L,NY,NX)
  micstt%CNO3B=trc_solcl_vr(ids_NO3B,L,NY,NX)
  micfor%RNH4EcoDmndSoilPrev =RNH4EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%RNO3EcoDmndSoilPrev =RNO3EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%RH2PO4EcoDmndSoilPrev =RH2PO4EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%RH1PO4EcoDmndSoilPrev =RH1PO4EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%VLWatMicP  =VLWatMicP(L,NY,NX)

  if(micfor%Lsurf)then
    micfor%SoilMicPMassLayer0=SoilMicPMassLayer(0,NY,NX)
  endif
  micfor%DiffusivitySolutEff(1:NPH)=DiffusivitySolutEff(1:NPH,L,NY,NX)
  micfor%FILM(1:NPH)=FILM(1:NPH,L,NY,NX)
  micfor%THETPM(1:NPH)=THETPM(1:NPH,L,NY,NX)
  micfor%VLWatMicPM(1:NPH)=VLWatMicPM(1:NPH,L,NY,NX)
  micfor%TortMicPM(1:NPH)=TortMicPM(1:NPH,L,NY,NX)
  micfor%VLsoiAirPM(1:NPH)=VLsoiAirPM(1:NPH,L,NY,NX)
  micfor%VLsoiAirP=VLsoiAirP(L,NY,NX)
  micstt%EPOC=EPOC(L,NY,NX)
  micstt%EHUM=EHUM(L,NY,NX)
  micstt%ZNH4B=trc_solml_vr(ids_NH4B,L,NY,NX)
  micstt%ZNH4S=trc_solml_vr(ids_NH4,L,NY,NX)
  micstt%ZNO3B=trc_solml_vr(ids_NO3B,L,NY,NX)
  micstt%ZNO3S=trc_solml_vr(ids_NO3,L,NY,NX)
  micstt%H1POB=trc_solml_vr(ids_H1PO4B,L,NY,NX)
  micstt%H1PO4=trc_solml_vr(ids_H1PO4,L,NY,NX)
  micstt%ZNO2B=trc_solml_vr(ids_NO2B,L,NY,NX)
  micstt%ZNO2S=trc_solml_vr(ids_NO2,L,NY,NX)
  micstt%H2POB=trc_solml_vr(ids_H2PO4B,L,NY,NX)
  micstt%H2PO4=trc_solml_vr(ids_H2PO4,L,NY,NX)
  micstt%CCO2S=trc_solcl_vr(idg_CO2,L,NY,NX)
  micstt%CNO2S=trc_solcl_vr(ids_NO2,L,NY,NX)
  micstt%CNO2B=trc_solcl_vr(ids_NO2B,L,NY,NX)
  micstt%CZ2OS=trc_solcl_vr(idg_N2O,L,NY,NX)
  micstt%Z2OS=trc_solml_vr(idg_N2O,L,NY,NX)
  micstt%COXYS=trc_solcl_vr(idg_O2,L,NY,NX)
  micstt%OXYS=trc_solml_vr(idg_O2,L,NY,NX)
  micstt%O2GSolubility=GasSolbility_vr(idg_O2,L,NY,NX)
  micstt%COXYG=trc_gascl_vr(idg_O2,L,NY,NX)
  micstt%CZ2GS=trc_solcl_vr(idg_N2,L,NY,NX)
  micstt%CH2GS=trc_solcl_vr(idg_H2,L,NY,NX)
  micstt%H2GS=trc_solml_vr(idg_H2,L,NY,NX)
  micstt%CCH4G=trc_gascl_vr(idg_CH4,L,NY,NX)
  micstt%CH4S=trc_solml_vr(idg_CH4,L,NY,NX)
  micstt%SCH4L=GasSolbility_vr(idg_CH4,L,NY,NX)
  micstt%ZNFN0=ZNFN0(L,NY,NX)
  micstt%ZNFNI=ZNFNI(L,NY,NX)
  micstt%FracBulkSOM(1:jcplx)=FracBulkSOM_vr(1:jcplx,L,NY,NX)
  micstt%DOM(idom_beg:idom_end,1:jcplx)=DOM(idom_beg:idom_end,1:jcplx,L,NY,NX)
  micstt%SorbedOM(idom_beg:idom_end,1:jcplx)=SorbedOM_vr(idom_beg:idom_end,1:jcplx,L,NY,NX)

  micstt%SolidOMAct(1:jsken,1:jcplx)=SolidOMAct_vr(1:jsken,1:jcplx,L,NY,NX)
  micstt%SolidOM(1:NumPlantChemElms,1:jsken,1:jcplx)=SolidOM_vr(1:NumPlantChemElms,1:jsken,1:jcplx,L,NY,NX)
  micstt%OMBioResdu(1:NumPlantChemElms,1:ndbiomcp,1:jcplx)=OMBioResdu_vr(1:NumPlantChemElms,1:ndbiomcp,1:jcplx,L,NY,NX)
  micstt%CNOSC(1:jsken,1:jcplx)=CNOSC(1:jsken,1:jcplx,L,NY,NX)
  micstt%CPOSC(1:jsken,1:jcplx)=CPOSC(1:jsken,1:jcplx,L,NY,NX)
  micstt%OMEheter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx)=OMEheter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,L,NY,NX)
  micstt%OMEauto(1:NumPlantChemElms,1:NumLiveAutoBioms)=OMEauto(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)
  if(.not.micfor%litrm)then
    micfor%AEC=AEC(L,NY,NX)
    micstt%OXYG=trc_gasml_vr(idg_O2,L,NY,NX)
  endif
  micflx%RNO2DmndSoilChemo=RNO2DmndSoilChemo_vr(L,NY,NX)
  micflx%RNO2DmndBandChemo=RNO2DmndBandChemo_vr(L,NY,NX)
  micflx%RNH4UptkSoilAutor(1:NumMicrobAutrophCmplx)=RNH4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNH4UptkBandAutor(1:NumMicrobAutrophCmplx)=RNH4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO3UptkSoilAutor(1:NumMicrobAutrophCmplx)=RNO3UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO3UptkBandAutor(1:NumMicrobAutrophCmplx)=RNO3UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH2PO4UptkSoilAutor(1:NumMicrobAutrophCmplx)=RH2PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH2PO4UptkBandAutor(1:NumMicrobAutrophCmplx)=RH2PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH1PO4UptkSoilAutor(1:NumMicrobAutrophCmplx)=RH1PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH1PO4UptkBandAutor(1:NumMicrobAutrophCmplx)=RH1PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RO2DmndAutort(1:NumMicrobAutrophCmplx)=RO2DmndAutort(1:NumMicrobAutrophCmplx,L,NY,NX)

  micflx%RNH4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)=RNH4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RNH4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)=RNH4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RNO3DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)=RNO3DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RNO3DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)=RNO3DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH2PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)=RH2PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH2PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)=RH2PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH1PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)=RH1PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH1PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)=RH1PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx)=RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)

  end subroutine MicAPISend

!------------------------------------------------------------------------------------------


  subroutine MicAPIRecv(L,NY,NX,litrM,micstt,micflx)
  implicit none
  integer, intent(in) :: L,NY,NX
  logical, intent(in) :: litrM
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(in) :: micflx
  integer :: NumMicbFunGroups, jcplx, NumMicrobAutrophCmplx
  NumMicrobAutrophCmplx = micpar%NumMicrobAutrophCmplx
  NumMicbFunGroups=micpar%NumMicbFunGroups
  jcplx=micpar%jcplx

  trcg_RMicbTransf_vr(idg_CO2,L,NY,NX) =micflx%RCO2NetUptkMicb
  trcg_RMicbTransf_vr(idg_CH4,L,NY,NX) =micflx%RCH4UptkAutor
  trcg_RMicbTransf_vr(idg_H2,L,NY,NX) =micflx%RH2NetUptkMicb
  trcg_RMicbTransf_vr(idg_O2,L,NY,NX)=micflx%RO2UptkMicb
  trcg_RMicbTransf_vr(idg_N2,L,NY,NX)  =micflx%RN2NetUptkMicb
  trcg_RMicbTransf_vr(idg_N2O,L,NY,NX)  =micflx%RN2ONetUptkMicb
  RNutMicbTransf_vr(ids_NH4,L,NY,NX) =micflx%RNH4MicbTransf_vr
  RNutMicbTransf_vr(ids_NO3,L,NY,NX) =micflx%RNO3MicbTransf_vr
  RNutMicbTransf_vr(ids_NO2,L,NY,NX) =micflx%RNO2MicbTransf_vr
  RNutMicbTransf_vr(ids_H2PO4,L,NY,NX) =micflx%RH2PO4MicbTransf_vr
  RNutMicbTransf_vr(ids_H1PO4,L,NY,NX) =micflx%RH1PO4MicbTransf_vr
  RNutMicbTransf_vr(ids_NH4B,L,NY,NX) =micflx%XNH4B
  RNutMicbTransf_vr(ids_NO3B,L,NY,NX) =micflx%XNO3B
  RNutMicbTransf_vr(ids_NO2B,L,NY,NX) =micflx%XNO2B
  RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX) =micflx%XH2BS
  RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX) =micflx%XH1BS
  Micb_N2Fixation_vr(L,NY,NX) =micflx%XN2GS
  RNO2DmndSoilChemo_vr(L,NY,NX)=micflx%RNO2DmndSoilChemo
  RNO2DmndBandChemo_vr(L,NY,NX)=micflx%RNO2DmndBandChemo
  NetNH4Mineralize_col(NY,NX)=NetNH4Mineralize_col(NY,NX)+micflx%NetNH4Mineralize
  NetPO4Mineralize_col(NY,NX)=NetPO4Mineralize_col(NY,NX)+micflx%NetPO4Mineralize_col
  RDOM_micb_flx(idom_doc,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_doc,1:jcplx)
  RDOM_micb_flx(idom_don,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_don,1:jcplx)
  RDOM_micb_flx(idom_dop,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_dop,1:jcplx)
  RDOM_micb_flx(idom_acetate,1:jcplx,L,NY,NX)=micflx%RDOM_micb_flx(idom_acetate,1:jcplx)

  RO2DmndAutort(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RO2DmndAutort(1:NumMicrobAutrophCmplx)
  RNH3OxidAutor(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNH3OxidAutor(1:NumMicrobAutrophCmplx)
  RNH3OxidAutorBand(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNH3OxidAutorBand(1:NumMicrobAutrophCmplx)
  RNO2OxidAutor(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNO2OxidAutor(1:NumMicrobAutrophCmplx)
  RNO2OxidAutorBand(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNO2OxidAutorBand(1:NumMicrobAutrophCmplx)
  RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx)
  RDOCUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RDOCUptkHeter(1:NumHetetrMicCmplx,1:jcplx)
  RAcetateUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RAcetateUptkHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3ReduxDmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNO3ReduxDmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3ReduxDmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNO3ReduxDmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO2DmndReduxSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNO2DmndReduxSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO2DmndReduxBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNO2DmndReduxBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RN2ODmndReduxHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RN2ODmndReduxHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNH4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNH4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNH4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNH4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNO3DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RNO3DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH2PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RH2PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH2PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RH2PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH1PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RH1PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH1PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)=micflx%RH1PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)

  SolidOM_vr(1:NumPlantChemElms,1:jsken,1:jcplx,L,NY,NX)=micstt%SolidOM(1:NumPlantChemElms,1:jsken,1:jcplx)
  SolidOMAct_vr(1:jsken,1:jcplx,L,NY,NX)=micstt%SolidOMAct(1:jsken,1:jcplx)

  if(litrm)then
    RNH4DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX)=micflx%RNH4DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RNO3DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX)=micflx%RNO3DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RH2PO4DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX)=micflx%RH2PO4DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RH1PO4DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX)=micflx%RH1PO4DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RNH4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)=micflx%RNH4UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RNO3UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)=micflx%RNO3UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RH2PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)=micflx%RH2PO4UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RH1PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)=micflx%RH1PO4UptkLitrAutor(1:NumMicrobAutrophCmplx)
    SolidOM_vr(ielmc,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)=micstt%OSC13U
    SolidOM_vr(ielmc,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSC14U
    SolidOM_vr(ielmc,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSC24U

    SolidOM_vr(ielmn,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)=micstt%OSN13U
    SolidOM_vr(ielmn,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSN14U
    SolidOM_vr(ielmn,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSN24U

    SolidOM_vr(ielmp,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)=micstt%OSP13U
    SolidOM_vr(ielmp,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSP14U
    SolidOM_vr(ielmp,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)=micstt%OSP24U
  endif

  RNH4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNH4UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RNH4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNH4UptkBandAutor(1:NumMicrobAutrophCmplx)
  RNO3UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNO3UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RNO3UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RNO3UptkBandAutor(1:NumMicrobAutrophCmplx)
  RH2PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RH2PO4UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RH2PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RH2PO4UptkBandAutor(1:NumMicrobAutrophCmplx)
  RH1PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RH1PO4UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RH1PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)=micflx%RH1PO4UptkBandAutor(1:NumMicrobAutrophCmplx)
  RO2UptkSoilM_vr(1:NPH,L,NY,NX)=micflx%RO2UptkSoilM(1:NPH)

  TSens4MicbGrwoth_vr(L,NY,NX)=micstt%TSens4MicbGrwoth
  VOLQ(L,NY,NX)=micstt%VOLQ
  TOQCK(L,NY,NX)=micstt%TOQCK
  ZNFNI(L,NY,NX)=micstt%ZNFNI
  FracBulkSOM_vr(1:jcplx,L,NY,NX)=micstt%FracBulkSOM(1:jcplx)
  DOM(idom_beg:idom_end,1:jcplx,L,NY,NX)=micstt%DOM(idom_beg:idom_end,1:jcplx)
  SorbedOM_vr(idom_beg:idom_end,1:jcplx,L,NY,NX)=micstt%SorbedOM(idom_beg:idom_end,1:jcplx)

  OMBioResdu_vr(1:NumPlantChemElms,1:ndbiomcp,1:jcplx,L,NY,NX)=micstt%OMBioResdu(1:NumPlantChemElms,1:ndbiomcp,1:jcplx)
  OMEheter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,L,NY,NX)=micstt%OMEheter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx)
  OMEauto(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)=micstt%OMEauto(1:NumPlantChemElms,1:NumLiveAutoBioms)

  end subroutine MicAPIRecv
end module MicBGCAPI
