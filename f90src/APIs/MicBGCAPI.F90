module MicBGCAPI

  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use SoilBGCNLayMod,       only: DownwardMixOM,          sumMicBiomLayL
  use SoilDisturbMod,       only: SOMRemovalByDisturbance
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicrobeDiagTypes,     only: Cumlate_Flux_Diag_type, Microbe_Diag_type
  use MicForcTypeMod,       only: micforctype
  use minimathmod,          only: AZMAX1,safe_adb,AZERO,AZERO1
  use EcoSiMParDataMod,     only: micpar
  use MicBGCMod,            only: SoilBGCOneLayer
  use EcosimConst,          only: LtHeatIceMelt,Tref
  use abortutils,           only: endrun
  use NumericalAuxMod
  use EcoSIMSolverPar  
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
implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  type(micforctype) :: micfor
  type(micsttype) :: micstt
  type(micfluxtype) :: micflx
  type(Microbe_Diag_type) :: nmicdiag

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
  call nmicdiag%Init()

  end subroutine MicAPI_Init
!------------------------------------------------------------------------------------------
  subroutine MicAPI_cleanup()

  implicit none


  call micfor%destroy()
  call micflx%destroy()
  call micstt%destroy()
  call nmicdiag%Destroy()

  end subroutine MicAPI_cleanup
!------------------------------------------------------------------------------------------
  subroutine MicrobeModel(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOIL BIOLOGICAL TRANSFORMATIONS
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: OrGM_beg(1:NumPlantChemElms)
  real(r8) :: dOrGM(1:NumPlantChemElms)
  real(r8) :: tdOrGM(1:NumPlantChemElms)
  integer :: L,NX,NY

!   begin_execution
  curI=I; curJ=J
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
!
!       VOLWZ=water volume used to calculate aqueous microbial
!       concentrations that drive microbial density effects on
!       decomposition
      D998: DO L=0,NL_col(NY,NX)
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          IF(L.EQ.0 .OR. L.GE.NU_col(NY,NX))THEN
             call sumMicBiomLayL(L,NY,NX,OrGM_beg)
             call MicBGC1Layer(I,J,L,NY,NX)
             call sumMicBiomLayL(L,NY,NX,dOrGM)
             dOrGM  = dOrGM-OrGM_beg
             tdOrGM = tdOrGM+dOrGM
          ELSE
            trcs_RMicbUptake_vr(idg_beg:idg_NH3-1,L,NY,NX)     = 0.0_r8
            RNut_MicbRelease_vr(ids_NH4B:ids_nuts_end,L,NY,NX) = 0.0_r8
            Micb_N2Fixation_vr(L,NY,NX)                        = 0.0_r8
          ENDIF
  
        ELSE
          trcs_RMicbUptake_vr(idg_beg:idg_NH3-1,L,NY,NX)     = 0.0_r8
          RNut_MicbRelease_vr(ids_NH4B:ids_nuts_end,L,NY,NX) = 0.0_r8
          Micb_N2Fixation_vr(L,NY,NX)                        = 0.0_r8
        ENDIF
      ENDDO D998
      
      DO L=0,NL_col(NY,NX)
        !     MIX LITTER C BETWEEN ADJACENT SOIL LAYERS L AND LL
        call DownwardMixOM(I,J,L,NY,NX,FracLitrMix_vr(L,NY,NX))
      ENDDO
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
  type(Cumlate_Flux_Diag_type) :: naqfdiag

  micfor%L=L
  call MicAPISend(I,J,L,NY,NX,micfor,micstt,micflx)
  
  call SoilBGCOneLayer(I,J,micfor,micstt,micflx,naqfdiag,nmicdiag)

  call MicAPIRecv(I,J,L,NY,NX,micfor%litrm,micstt,micflx,naqfdiag,nmicdiag)
  
  end subroutine Micbgc1Layer
!------------------------------------------------------------------------------------------

  subroutine MicAPISend(I,J,L,NY,NX,micfor,micstt,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: L,NY,NX
  type(micforctype), intent(inout) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx

  integer :: NumMicbFunGrupsPerCmplx, jcplx, k_POM, k_humus, idom, K,KL
  integer :: ndbiomcp, nlbiomcp, NumMicrobAutrophCmplx, NumHetetr1MicCmplx,NE
  
  NumMicbFunGrupsPerCmplx = micpar%NumMicbFunGrupsPerCmplx
  jcplx                   = micpar%jcplx

  ndbiomcp              = micpar%ndbiomcp
  nlbiomcp              = micpar%nlbiomcp
  NumMicrobAutrophCmplx = micpar%NumMicrobAutrophCmplx
  NumHetetr1MicCmplx     = micpar%NumHetetr1MicCmplx
  k_humus               = micpar%k_humus
  k_POM                 = micpar%k_POM

  !is it litter layer?
  micfor%litrm=(L==0)
  !is it surface layer
  micfor%Lsurf=(L==NU_col(NY,NX))

  micfor%VOLW0               = VLWatMicP_vr(0,NY,NX)
  micfor%ZERO                = ZERO
  micfor%CCH4E               = AtmGasCgperm3_col(idg_CH4,NY,NX)
  micfor%COXYE               = AtmGasCgperm3_col(idg_O2,NY,NX)
  micfor%O2_irrig_conc       = trcg_irrig_mole_conc_col(idg_O2,NY,NX)
  micfor%O2_rain_conc        = trcg_rain_mole_conc_col(idg_O2,NY,NX)
  micfor%Irrig2LitRSurf_col  = Irrig2LitRSurf_col(NY,NX)
  micfor%Rain2LitRSurf       = Rain2LitRSurf_col(NY,NX)
  micfor%TempOffset          = TempOffset_col(NY,NX)
  micfor%VLitR               = VLitR_col(NY,NX)
  micfor%VWatLitRHoldCapcity = VWatLitRHoldCapcity_col(NY,NX)
  micfor%ZEROS2              = ZEROS2(NY,NX)
  micfor%ZEROS               = ZEROS(NY,NX)
  micfor%VLSoilMicP          = VLSoilMicP_vr(L,NY,NX)
  micfor%THETY               = SoilWatAirDry_vr(L,NY,NX)
  micfor%POROS               = POROS_vr(L,NY,NX)
  micfor%FieldCapacity       = FieldCapacity_vr(L,NY,NX)

  micfor%THETW             = THETW_vr(L,NY,NX)
  micfor%PH                = PH_vr(L,NY,NX)
  micfor%SoilMicPMassLayer = VLSoilMicPMass_vr(L,NY,NX)
  micfor%VLSoilPoreMicP    = VLSoilPoreMicP_vr(L,NY,NX)
  micfor%TScal4Difsvity    = TScal4Difsvity_vr(L,NY,NX)
  micfor%VLNOB             = trcs_VLN_vr(ids_NO3B,L,NY,NX)
  micfor%VLNO3             = trcs_VLN_vr(ids_NO3,L,NY,NX)
  micfor%VLNH4             = trcs_VLN_vr(ids_NH4,L,NY,NX)
  micfor%VLNHB             = trcs_VLN_vr(ids_NH4B,L,NY,NX)
  micfor%VLPO4             = trcs_VLN_vr(ids_H1PO4,L,NY,NX)
  micfor%VLPOB             = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

  if(micfor%litrm)then
    KL=micpar%NumOfLitrCmplxs  
    micfor%TKS            = FracSurfByLitR_col(NY,NX)*TKS_vr(0,NY,NX)+(1._r8-FracSurfByLitR_col(NY,NX))*TKS_vr(NU_col(NY,NX),NY,NX)
    micfor%PSISoilMatricP = FracSurfByLitR_col(NY,NX)*PSISoilMatricP_vr(0,NY,NX)+ &
      (1._r8-FracSurfByLitR_col(NY,NX))*PSISoilMatricP_vr(NU_col(NY,NX),NY,NX)
  else
    KL=jcplx
    micfor%TKS            = TKS_vr(L,NY,NX)
    micfor%PSISoilMatricP = PSISoilMatricP_vr(L,NY,NX)
  endif  

  if (micfor%TKS<Tref)then
    micfor%PSISoilMatricP  = LtHeatIceMelt*(micfor%TKS-Tref)/micfor%TKS
  endif
  
  micfor%O2AquaDiffusvity         = SoluteDifusvty_vr(idg_O2,L,NY,NX)
  micfor%ORGC                     = SoilOrgM_vr(ielmc,L,NY,NX)
  micfor%RNO2EcoUptkSoilPrev      = RNO2EcoUptkSoilPrev_vr(L,NY,NX)
  micfor%RN2OEcoUptkSoilPrev      = RN2OEcoUptkSoilPrev_vr(L,NY,NX)
  micfor%RNO2EcoUptkBandPrev      = RNO2EcoUptkBandPrev_vr(L,NY,NX)
  micfor%RO2EcoDmndPrev           = RO2EcoDmndPrev_vr(L,NY,NX)
  micfor%RO2GasXchangePrev        = RGasTranspFlxPrev_vr(idg_O2,L,NY,NX)
  micfor%RCH4GasXchangePrev       = RGasTranspFlxPrev_vr(idg_CH4,L,NY,NX)
  micfor%RCH4PhysexchPrev         = RCH4PhysexchPrev_vr(L,NY,NX)
  micfor%RNH4EcoDmndBandPrev      = RNH4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RNO3EcoDmndBandPrev      = RNO3EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RH2PO4EcoDmndBandPrev    = RH2PO4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RH1PO4EcoDmndBandPrev    = RH1PO4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RO2AquaXchangePrev       = RO2AquaSourcePrev_vr(L,NY,NX)
  micfor%RDOMEcoDmndPrev(1:KL) = RDOMEcoDmndPrev_vr(1:KL,L,NY,NX)

!  write(201,*)I+J/24.,L,RDOMEcoDmndPrev_vr(1:KL,L,NY,NX)
  micfor%RAcetateEcoDmndPrev(1:KL) = RAcetateEcoDmndPrev_vr(1:KL,L,NY,NX)

  !litter layer is modeled
  if(micfor%litrm)then
    micstt%ZNH4TU = AZMAX1(trcs_solml_vr(ids_NH4,NU_col(NY,NX),NY,NX))+AZMAX1(trcs_solml_vr(ids_NH4B,NU_col(NY,NX),NY,NX))
    micstt%ZNO3TU = AZMAX1(trcs_solml_vr(ids_NO3,NU_col(NY,NX),NY,NX))+AZMAX1(trcs_solml_vr(ids_NO3B,NU_col(NY,NX),NY,NX))
    micstt%H1P4TU = AZMAX1(trcs_solml_vr(ids_H1PO4,NU_col(NY,NX),NY,NX))+AZMAX1(trcs_solml_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX))
    micstt%H2P4TU = AZMAX1(trcs_solml_vr(ids_H2PO4,NU_col(NY,NX),NY,NX))+AZMAX1(trcs_solml_vr(ids_H2PO4B,NU_col(NY,NX),NY,NX))

    micstt%CNH4BU  = trc_solcl_vr(ids_NH4B,NU_col(NY,NX),NY,NX)
    micstt%CNH4SU  = trc_solcl_vr(ids_NH4,NU_col(NY,NX),NY,NX)
    micstt%CH2P4U  = trc_solcl_vr(ids_H2PO4,NU_col(NY,NX),NY,NX)
    micstt%CH2P4BU = trc_solcl_vr(ids_H2PO4B,NU_col(NY,NX),NY,NX)
    micstt%CH1P4U  = trc_solcl_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
    micstt%CH1P4BU = trc_solcl_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
    micstt%CNO3SU  = trc_solcl_vr(ids_NO3,NU_col(NY,NX),NY,NX)
    micstt%CNO3BU  = trc_solcl_vr(ids_NO3B,NU_col(NY,NX),NY,NX)

    DO NE=1,NumPlantChemElms
      micstt%SOMPomProtein(NE)  = SolidOM_vr(NE,micpar%iprotein,micpar%k_POM,NU_col(NY,NX),NY,NX)
      micstt%SOMHumProtein(NE)  = SolidOM_vr(NE,micpar%iprotein,micpar%k_humus,NU_col(NY,NX),NY,NX)
      micstt%SOMHumCarbohyd(NE) = SolidOM_vr(NE,micpar%icarbhyro,micpar%k_humus,NU_col(NY,NX),NY,NX)
    ENDDO
    micfor%RNH4EcoDmndLitrPrev       = RNH4EcoDmndSoilPrev_vr(NU_col(NY,NX),NY,NX)
    micfor%RNO3EcoDmndLitrPrev       = RNO3EcoDmndSoilPrev_vr(NU_col(NY,NX),NY,NX)
    micfor%RH2PO4EcoDmndLitrPrev     = RH2PO4EcoDmndSoilPrev_vr(NU_col(NY,NX),NY,NX)
    micfor%RH1PO4EcoDmndLitrPrev     = RH1PO4EcoDmndSoilPrev_vr(NU_col(NY,NX),NY,NX)
    micfor%VOLWU                     = VLWatMicP_vr(NU_col(NY,NX),NY,NX)
    micfor%ElmAllocmatMicrblitr2POMU = ElmAllocmatMicrblitr2POM_vr(1:2,NU_col(NY,NX),NY,NX)

    micflx%RNH4DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:KL)   = RNH4DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX)
    micflx%RNO3DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:KL)   = RNO3DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX)
    micflx%RH2PO4DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:KL) = RH2PO4DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX)
    micflx%RH1PO4DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:KL) = RH1PO4DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX)
    micflx%RNH4UptkLitrAutorPrev(1:NumMicrobAutrophCmplx)     = RNH4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)
    micflx%RNO3UptkLitrAutorPrev(1:NumMicrobAutrophCmplx)     = RNO3UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)
    micflx%RH2PO4UptkLitrAutorPrev(1:NumMicrobAutrophCmplx)   = RH2PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)
    micflx%RH1PO4UptkLitrAutorPrev(1:NumMicrobAutrophCmplx)   = RH1PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)

  else
    micfor%ElmAllocmatMicrblitr2POM =ElmAllocmatMicrblitr2POM_vr(1:2,L,NY,NX)
  endif

  micstt%CNH4B                 = trc_solcl_vr(ids_NH4B,L,NY,NX)
  micstt%CNH4S                 = trc_solcl_vr(ids_NH4,L,NY,NX)
  micstt%CH2P4                 = trc_solcl_vr(ids_H2PO4,L,NY,NX)
  micstt%CH2P4B                = trc_solcl_vr(ids_H2PO4B,L,NY,NX)
  micstt%CH1P4                 = trc_solcl_vr(ids_H1PO4,L,NY,NX)
  micstt%CH1P4B                = trc_solcl_vr(ids_H1PO4B,L,NY,NX)
  micstt%CNO3S                 = trc_solcl_vr(ids_NO3,L,NY,NX)
  micstt%CNO3B                 = trc_solcl_vr(ids_NO3B,L,NY,NX)
  micfor%RNH4EcoDmndSoilPrev   = RNH4EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%RNO3EcoDmndSoilPrev   = RNO3EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%RH2PO4EcoDmndSoilPrev = RH2PO4EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%RH1PO4EcoDmndSoilPrev = RH1PO4EcoDmndSoilPrev_vr(L,NY,NX)
  micfor%VLWatMicP             = VLWatMicP_vr(L,NY,NX)

  if(micfor%Lsurf)then
    micfor%SoilMicPMassLayer0=VLSoilMicPMass_vr(0,NY,NX)
  endif
  micfor%DiffusivitySolutEff(1:NPH) = DiffusivitySolutEffM_vr(1:NPH,L,NY,NX)
  micfor%FILM(1:NPH)                = FILMM_vr(1:NPH,L,NY,NX)
  micfor%THETPM(1:NPH)              = FracAirFilledSoilPoreM_vr(1:NPH,L,NY,NX)
  micfor%VLWatMicPM(1:NPH)          = VLWatMicPM_vr(1:NPH,L,NY,NX)
  micfor%TortMicPM(1:NPH)           = TortMicPM_vr(1:NPH,L,NY,NX)
  micfor%VLsoiAirPM(1:NPH)          = VLsoiAirPM_vr(1:NPH,L,NY,NX)
  micfor%VLsoiAirP                  = VLsoiAirP_vr(L,NY,NX)
  micstt%EPOC                       = EPOC_vr(L,NY,NX)
  micstt%EHUM                       = EHUM_vr(L,NY,NX)
  micstt%ZNH4B             = AZMAX1(trcs_solml_vr(ids_NH4B,L,NY,NX)-trcs_solml_drib_vr(ids_NH4B,L,NY,NX))
  micstt%ZNH4S             = AZMAX1(trcs_solml_vr(ids_NH4,L,NY,NX)-trcs_solml_drib_vr(ids_NH4,L,NY,NX))
  micstt%ZNO3B             = AZMAX1(trcs_solml_vr(ids_NO3B,L,NY,NX)-trcs_solml_drib_vr(ids_NO3B,L,NY,NX))
  micstt%ZNO3S             = AZMAX1(trcs_solml_vr(ids_NO3,L,NY,NX)-trcs_solml_drib_vr(ids_NO3,L,NY,NX))
  micstt%H1POB             = AZMAX1(trcs_solml_vr(ids_H1PO4B,L,NY,NX)-trcs_solml_drib_vr(ids_H1PO4B,L,NY,NX))
  micstt%H1PO4             = AZMAX1(trcs_solml_vr(ids_H1PO4,L,NY,NX)-trcs_solml_drib_vr(ids_H1PO4,L,NY,NX))
  micstt%ZNO2B             = AZMAX1(trcs_solml_vr(ids_NO2B,L,NY,NX)-trcs_solml_drib_vr(ids_NO2B,L,NY,NX))
  micstt%ZNO2S             = AZMAX1(trcs_solml_vr(ids_NO2,L,NY,NX)-trcs_solml_drib_vr(ids_NO2,L,NY,NX))
  micstt%H2POB             = AZMAX1(trcs_solml_vr(ids_H2PO4B,L,NY,NX)-trcs_solml_drib_vr(ids_H2PO4B,L,NY,NX))
  micstt%H2PO4             = AZMAX1(trcs_solml_vr(ids_H2PO4,L,NY,NX)-trcs_solml_drib_vr(ids_H2PO4,L,NY,NX))
  micstt%Z2OS              = AZMAX1(trcs_solml_vr(idg_N2O,L,NY,NX)-trcs_solml_drib_vr(idg_N2O,L,NY,NX))  
  micstt%CH4S              = AZMAX1(trcs_solml_vr(idg_CH4,L,NY,NX)-trcs_solml_drib_vr(idg_CH4,L,NY,NX))  
  micstt%OXYS              = AZMAX1(trcs_solml_vr(idg_O2,L,NY,NX)-trcs_solml_drib_vr(idg_O2,L,NY,NX))    
  micstt%H2GS              = AZMAX1(trcs_solml_vr(idg_H2,L,NY,NX)-trcs_solml_drib_vr(idg_H2,L,NY,NX))  
  micstt%CCO2S             = AZMAX1(trc_solcl_vr(idg_CO2,L,NY,NX))  
  micstt%CNO2S             = AZMAX1(trc_solcl_vr(ids_NO2,L,NY,NX))  
  micstt%CNO2B             = AZMAX1(trc_solcl_vr(ids_NO2B,L,NY,NX))  
  micstt%CZ2OS             = AZMAX1(trc_solcl_vr(idg_N2O,L,NY,NX))  
  micstt%COXYS             = AZMAX1(trc_solcl_vr(idg_O2,L,NY,NX))  
  micstt%COXYG             = AZMAX1(trcg_gascl_vr(idg_O2,L,NY,NX))  
  micstt%CZ2GS             = AZMAX1(trc_solcl_vr(idg_N2,L,NY,NX))
  micstt%CH2GS             = AZMAX1(trc_solcl_vr(idg_H2,L,NY,NX))
  micstt%CCH4G             = AZMAX1(trcg_gascl_vr(idg_CH4,L,NY,NX))
  micstt%Lay               = L

  FermOXYI_vr(L,NY,NX)  = 1.0_r8-1.0_r8/(1.0_r8+EXP(1.0_r8*AMAX1(-micstt%COXYS+2.5_r8,-50._r8)))

!  write(115,*)I+J/24.,L,micstt%COXYG,micstt%COXYS,VLsoiAirPM(1,L,NY,NX)
  micstt%O2GSolubility     = GasSolbility_vr(idg_O2,L,NY,NX)  
  micstt%CH4AquaSolubility = GasSolbility_vr(idg_CH4,L,NY,NX)
  micstt%ZNFN0             = ZNFN0_vr(L,NY,NX)
  micstt%ZNFNI             = ZNFNI_vr(L,NY,NX)

  if(.not.micfor%litrm)then
    micfor%AEC  = AEC_vr(L,NY,NX)
    micstt%OXYG = AZMAX1(trcg_gasml_vr(idg_O2,L,NY,NX))
  ENDIF 

  DO K=1,KL
    DO idom=idom_beg,idom_end  
      micstt%DOM(idom,K)           = AZMAX1(DOM_MicP_vr(idom,K,L,NY,NX))
      micstt%DOM_MicP_drib(idom,K) = DOM_MicP_drib_vr(idom,K,L,NY,NX)
    ENDDO
  ENDDO

  micstt%SorbedOM(idom_beg:idom_end,1:KL)               = SorbedOM_vr(idom_beg:idom_end,1:KL,L,NY,NX)
  micstt%SolidOMAct(1:jsken,1:KL)                       = SolidOMAct_vr(1:jsken,1:KL,L,NY,NX)
  micstt%SolidOM(1:NumPlantChemElms,1:jsken,1:KL)       = SolidOM_vr(1:NumPlantChemElms,1:jsken,1:KL,L,NY,NX)
  micstt%OMBioResdu(1:NumPlantChemElms,1:ndbiomcp,1:KL) = OMBioResdu_vr(1:NumPlantChemElms,1:ndbiomcp,1:KL,L,NY,NX)
  micstt%CNOSC(1:jsken,1:KL)                            = CNOSC_vr(1:jsken,1:KL,L,NY,NX)
  micstt%CPOSC(1:jsken,1:KL)                            = CPOSC_vr(1:jsken,1:KL,L,NY,NX)
  micstt%mBiomeHeter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:KL)=mBiomeHeter_vr(1:NumPlantChemElms,1:NumLiveHeterBioms,1:KL,L,NY,NX)

  micstt%mBiomeAutor(1:NumPlantChemElms,1:NumLiveAutoBioms)=mBiomeAutor_vr(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)

  micflx%RNO2DmndSoilChemoPrev=RNO2DmndSoilChemo_vr(L,NY,NX)
  micflx%RNO2DmndBandChemoPrev=RNO2DmndBandChemo_vr(L,NY,NX)
  micflx%RNH4UptkSoilAutorPrev(1:NumMicrobAutrophCmplx)   = RNH4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNH4UptkBandAutorPrev(1:NumMicrobAutrophCmplx)   = RNH4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO3UptkSoilAutorPrev(1:NumMicrobAutrophCmplx)   = RNO3UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO3UptkBandAutorPrev(1:NumMicrobAutrophCmplx)   = RNO3UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH2PO4UptkSoilAutorPrev(1:NumMicrobAutrophCmplx) = RH2PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH2PO4UptkBandAutorPrev(1:NumMicrobAutrophCmplx) = RH2PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH1PO4UptkSoilAutorPrev(1:NumMicrobAutrophCmplx) = RH1PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH1PO4UptkBandAutorPrev(1:NumMicrobAutrophCmplx) = RH1PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RO2DmndAutortPrev(1:NumMicrobAutrophCmplx)       = RO2DmndAutort_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO2OxidAutorPrev(1:NumMicrobAutrophCmplx)       = RNO2OxidAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO2OxidAutorBandPrev(1:NumMicrobAutrophCmplx)   = RNO2OxidAutorBand_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNH3OxidAutorPrev(1:NumMicrobAutrophCmplx) =RNH3OxidAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)       
  micflx%RNH3OxidAutorBandPrev(1:NumMicrobAutrophCmplx)=RNH3OxidAutorBand_vr(1:NumMicrobAutrophCmplx,L,NY,NX)

  micflx%RNH4DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:KL)      = RNH4DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RNH4DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:KL)      = RNH4DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RNO3DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:KL)      = RNO3DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RNO3DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:KL)      = RNO3DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RH2PO4DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:KL)    = RH2PO4DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RH2PO4DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:KL)    = RH2PO4DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RH1PO4DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:KL)    = RH1PO4DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RH1PO4DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:KL)    = RH1PO4DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RO2DmndHetertPrev(1:NumHetetr1MicCmplx,1:KL)          = RO2DmndHetert_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RDOCUptkHeterPrev(1:NumHetetr1MicCmplx,1:KL)          = RDOCUptkHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RAcetateUptkHeterPrev(1:NumHetetr1MicCmplx,1:KL)      = RAcetateUptkHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RNO3ReduxDmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:KL) = RNO3ReduxDmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RNO3ReduxDmndBandHeterPrev(1:NumHetetr1MicCmplx,1:KL) = RNO3ReduxDmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RNO2DmndReduxSoilHeterPrev(1:NumHetetr1MicCmplx,1:KL) = RNO2DmndReduxSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RNO2DmndReduxBandHeterPrev(1:NumHetetr1MicCmplx,1:KL) = RNO2DmndReduxBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  micflx%RN2ODmndReduxHeterPrev(1:NumHetetr1MicCmplx,1:KL)     = RN2ODmndReduxHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)
  end subroutine MicAPISend

!------------------------------------------------------------------------------------------


  subroutine MicAPIRecv(I,J,L,NY,NX,litrM,micstt,micflx,naqfdiag,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: L,NY,NX
  logical, intent(in) :: litrM
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(in) :: micflx
  type(Cumlate_Flux_Diag_type), intent(in) :: naqfdiag
  type(Microbe_Diag_type), intent(in) :: nmicdiag

  integer :: NumMicbFunGrupsPerCmplx, jcplx, NumMicrobAutrophCmplx
  integer :: NE,idom,K,idg,KL,NN
  
  NumMicrobAutrophCmplx = micpar%NumMicrobAutrophCmplx
  NumMicbFunGrupsPerCmplx=micpar%NumMicbFunGrupsPerCmplx
  jcplx=micpar%jcplx
  if(litrM)then
    KL=micpar%NumOfLitrCmplxs
  else
    KL=jcplx
  endif
  RCH4ProdAcetcl_vr(L,NY,NX)                     = naqfdiag%tCH4ProdAceto
  RCH4ProdHydrog_vr(L,NY,NX)                     = naqfdiag%tCH4ProdH2
  RCH4Oxi_aero_vr(L,NY,NX)                       = naqfdiag%tCH4OxiAero
  RFerment_vr(L,NY,NX)                           = naqfdiag%tCResp4H2Prod
  RNH3oxi_vr(L,NY,NX)                            = naqfdiag%tRNH3Oxi
  RN2ODeniProd_vr(L,NY,NX)                       = naqfdiag%TDeniReduxNO2Soil+naqfdiag%TDeniReduxNO2Band
  RN2OChemoProd_vr(L,NY,NX)                      = naqfdiag%RN2OProdSoilChemo+naqfdiag%RN2OProdBandChemo
  RN2ONitProd_vr(L,NY,NX)                        = naqfdiag%TNitReduxNO2Soil+naqfdiag%TNitReduxNO2Band
  RN2ORedux_vr(L,NY,NX)                          = naqfdiag%TReduxN2O
  OxyDecompLimiter_vr(L,NY,NX)                   = safe_adb(naqfdiag%tRO2UptkHeterG,naqfdiag%tRO2DmndHeterG)
  RO2DecompUptk_vr(L,NY,NX)                      = naqfdiag%tRO2UptkHeterG
  tRHydlySOM_vr(1:NumPlantChemElms,L,NY,NX)      = micflx%tRHydlySOM
  tRHydlyBioReSOM_vr(1:NumPlantChemElms,L,NY,NX) = micflx%tRHydlyBioReSOM
  tRHydlySoprtOM_vr(1:NumPlantChemElms,L,NY,NX)  = micflx%tRHydlySoprtOM
  ROQC4HeterMicActCmpK_vr(1:KL,L,NY,NX)          = nmicdiag%ROQC4HeterMicActCmpK(1:KL)
  RHydrolysisScalCmpK_vr(1:KL,L,NY,NX) = nmicdiag%RHydrolysisScalCmpK(1:KL)
  trcs_RMicbUptake_vr(idg_CO2,L,NY,NX) = micflx%RCO2NetUptkMicb
  trcs_RMicbUptake_vr(idg_CH4,L,NY,NX) = naqfdiag%tCH4OxiAero-naqfdiag%tCH4ProdAceto-naqfdiag%tCH4ProdH2
  trcs_RMicbUptake_vr(idg_H2,L,NY,NX)  = micflx%RH2NetUptkMicb
  trcs_RMicbUptake_vr(idg_O2,L,NY,NX)  = micflx%RO2UptkMicb
  trcs_RMicbUptake_vr(idg_N2,L,NY,NX)  = micflx%RN2NetUptkMicb+micflx%MicrbN2Fix
  trcs_RMicbUptake_vr(idg_N2O,L,NY,NX) = micflx%RN2ONetUptkMicb
  tRespGrossHeterUlm_vr(L,NY,NX)       = naqfdiag%tRespGrossHeterUlm
  tRespGrossHeter_vr(L,NY,NX)          = naqfdiag%tRespGrossHeter
  !artificial source
  if(L>0 .and. L<5)then
    trcs_RMicbUptake_vr(idg_Ar,L,NY,NX)  = -0.01_r8
  else
    trcs_RMicbUptake_vr(idg_Ar,L,NY,NX)  = 0._r8
  endif  

  RNut_MicbRelease_vr(ids_NH4,L,NY,NX)    = micflx%RNH4MicbReliz2Soil
  RNut_MicbRelease_vr(ids_NO3,L,NY,NX)    = micflx%RNO3MicbReliz2Soil
  RNut_MicbRelease_vr(ids_NO2,L,NY,NX)    = micflx%RNO2MicbReliz2Soil
  RNut_MicbRelease_vr(ids_H2PO4,L,NY,NX)  = micflx%RH2PO4MicbReliz2Soil
  RNut_MicbRelease_vr(ids_H1PO4,L,NY,NX)  = micflx%RH1PO4MicbReliz2Soil
  RNut_MicbRelease_vr(ids_NH4B,L,NY,NX)   = micflx%RNH4MicbReliz2Band
  RNut_MicbRelease_vr(ids_NO3B,L,NY,NX)   = micflx%RNO3MicbReliz2Band
  RNut_MicbRelease_vr(ids_NO2B,L,NY,NX)   = micflx%RNO2MicbReliz2Band
  RNut_MicbRelease_vr(ids_H2PO4B,L,NY,NX) = micflx%RH2PO4MicbReliz2Band
  RNut_MicbRelease_vr(ids_H1PO4B,L,NY,NX) = micflx%RH1PO4MicbReliz2Band
  
  TempSensDecomp_vr(L,NY,NX)            = nmicdiag%TSensGrowth
  MoistSensDecomp_vr(L,NY,NX)           = AZMAX1(nmicdiag%WatStressMicb)
  Micb_N2Fixation_vr(L,NY,NX)           = micflx%MicrbN2Fix    
  RNO2DmndSoilChemo_vr(L,NY,NX)         = micflx%RNO2DmndSoilChemo
  RNO2DmndBandChemo_vr(L,NY,NX)         = micflx%RNO2DmndBandChemo
  NetNH4Mineralize_CumYr_col(NY,NX)     = NetNH4Mineralize_CumYr_col(NY,NX)+micflx%NetNH4Mineralize
  NetPO4Mineralize_CumYr_col(NY,NX)     = NetPO4Mineralize_CumYr_col(NY,NX)+micflx%NetPO4Mineralize

  DO idom=idom_beg,idom_end
    REcoDOMProd_vr(idom,1:KL,L,NY,NX)=micflx%REcoDOMProd(idom,1:KL)    
  ENDDO
  RDOMMicProd_vr(idom_beg:idom_end,1:KL,L,NY,NX)=REcoDOMProd_vr(idom_beg:idom_end,1:KL,L,NY,NX)
 
  RO2DmndAutort_vr(1:NumMicrobAutrophCmplx,L,NY,NX)            = micflx%RO2DmndAutort(1:NumMicrobAutrophCmplx)
  RNH3OxidAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)            = micflx%RNH3OxidAutor(1:NumMicrobAutrophCmplx)
  RNH3OxidAutorBand_vr(1:NumMicrobAutrophCmplx,L,NY,NX)        = micflx%RNH3OxidAutorBand(1:NumMicrobAutrophCmplx)
  RNO2OxidAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)            = micflx%RNO2OxidAutor(1:NumMicrobAutrophCmplx)
  RNO2OxidAutorBand_vr(1:NumMicrobAutrophCmplx,L,NY,NX)        = micflx%RNO2OxidAutorBand(1:NumMicrobAutrophCmplx)
  RO2DmndHetert_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)          = micflx%RO2DmndHetert(1:NumHetetr1MicCmplx,1:KL)
  RDOCUptkHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)          = micflx%RDOCUptkHeter(1:NumHetetr1MicCmplx,1:KL)
  RAcetateUptkHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)      = micflx%RAcetateUptkHeter(1:NumHetetr1MicCmplx,1:KL)
  RNO3ReduxDmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX) = micflx%RNO3ReduxDmndSoilHeter(1:NumHetetr1MicCmplx,1:KL)
  RNO3ReduxDmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX) = micflx%RNO3ReduxDmndBandHeter(1:NumHetetr1MicCmplx,1:KL)
  RNO2DmndReduxSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX) = micflx%RNO2DmndReduxSoilHeter(1:NumHetetr1MicCmplx,1:KL)
  RNO2DmndReduxBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX) = micflx%RNO2DmndReduxBandHeter(1:NumHetetr1MicCmplx,1:KL)
  RN2ODmndReduxHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)     = micflx%RN2ODmndReduxHeter(1:NumHetetr1MicCmplx,1:KL)
  RNH4DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)      = micflx%RNH4DmndSoilHeter(1:NumHetetr1MicCmplx,1:KL)
  RNH4DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)      = micflx%RNH4DmndBandHeter(1:NumHetetr1MicCmplx,1:KL)
  RNO3DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)      = micflx%RNO3DmndSoilHeter(1:NumHetetr1MicCmplx,1:KL)
  RNO3DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)      = micflx%RNO3DmndBandHeter(1:NumHetetr1MicCmplx,1:KL)
  RH2PO4DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)    = micflx%RH2PO4DmndSoilHeter(1:NumHetetr1MicCmplx,1:KL)
  RH2PO4DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)    = micflx%RH2PO4DmndBandHeter(1:NumHetetr1MicCmplx,1:KL)
  RH1PO4DmndSoilHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)    = micflx%RH1PO4DmndSoilHeter(1:NumHetetr1MicCmplx,1:KL)
  RH1PO4DmndBandHeter_vr(1:NumHetetr1MicCmplx,1:KL,L,NY,NX)    = micflx%RH1PO4DmndBandHeter(1:NumHetetr1MicCmplx,1:KL)
  
  SolidOM_vr(1:NumPlantChemElms,1:jsken,1:KL,L,NY,NX) = micstt%SolidOM(1:NumPlantChemElms,1:jsken,1:KL)

  SolidOMAct_vr(1:jsken,1:KL,L,NY,NX)                 = micstt%SolidOMAct(1:jsken,1:KL)
  TSolidOMActC_vr(L,NY,NX)                            = micstt%TSolidOMActC
  TSolidOMC_vr(L,NY,NX)                               = micstt%TSolidOMC
  tOMActC_vr(L,NY,NX)                                 = micstt%tOMActC

  if(litrM)then
    RNH4DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX)   = micflx%RNH4DmndLitrHeter(1:NumHetetr1MicCmplx,1:KL)
    RNO3DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX)   = micflx%RNO3DmndLitrHeter(1:NumHetetr1MicCmplx,1:KL)
    RH2PO4DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX) = micflx%RH2PO4DmndLitrHeter(1:NumHetetr1MicCmplx,1:KL)
    RH1PO4DmndLitrHeter_col(1:NumHetetr1MicCmplx,1:KL,NY,NX) = micflx%RH1PO4DmndLitrHeter(1:NumHetetr1MicCmplx,1:KL)
    RNH4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)       = micflx%RNH4UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RNO3UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)       = micflx%RNO3UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RH2PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)     = micflx%RH2PO4UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RH1PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)     = micflx%RH1PO4UptkLitrAutor(1:NumMicrobAutrophCmplx)

    DO NE=1,NumPlantChemElms
      SolidOM_vr(NE,micpar%iprotein,micpar%k_POM,NU_col(NY,NX),NY,NX)    = micstt%SOMPomProtein(NE)
      SolidOM_vr(NE,micpar%iprotein,micpar%k_humus,NU_col(NY,NX),NY,NX)  = micstt%SOMHumProtein(NE)
      SolidOM_vr(NE,micpar%icarbhyro,micpar%k_humus,NU_col(NY,NX),NY,NX) = micstt%SOMHumCarbohyd(NE)
    ENDDO
  endif

  RNH4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)   = micflx%RNH4UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RNH4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)   = micflx%RNH4UptkBandAutor(1:NumMicrobAutrophCmplx)
  RNO3UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)   = micflx%RNO3UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RNO3UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)   = micflx%RNO3UptkBandAutor(1:NumMicrobAutrophCmplx)
  RH2PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX) = micflx%RH2PO4UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RH2PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX) = micflx%RH2PO4UptkBandAutor(1:NumMicrobAutrophCmplx)
  RH1PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX) = micflx%RH1PO4UptkSoilAutor(1:NumMicrobAutrophCmplx)
  RH1PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX) = micflx%RH1PO4UptkBandAutor(1:NumMicrobAutrophCmplx)
  DO NN=1,NPH
    REcoUptkSoilO2M_vr(NN,L,NY,NX)                       = REcoUptkSoilO2M_vr(NN,L,NY,NX)+micflx%REcoUptkSoilO2M(NN)
  ENDDO

  call nmicdiag%Summary(micpar%mid_Aerob_HeteroBacter,AeroBact_PrimeS_lim_vr(L,NY,NX))
  call nmicdiag%Summary(micpar%mid_Aerob_Fungi,AeroFung_PrimeS_lim_vr(L,NY,NX))

  TSens4MicbGrwoth_vr(L,NY,NX)                   = micstt%TSens4MicbGrwoth
  VWatMicrobAct_vr(L,NY,NX)                      = micstt%VWatMicrobAct
  TMicHeterActivity_vr(L,NY,NX)                  = micstt%TMicHeterActivity
  ZNFNI_vr(L,NY,NX)                              = micstt%ZNFNI

  do K=1,KL
    FracBulkSOMC_vr(K,L,NY,NX)    = micstt%FracBulkSOMC(K)
    RHydlySOCK_vr(K,L,NY,NX)      = micflx%RHydlySOCK(K)
    DO idom=idom_beg,idom_end
      DOM_MicP_vr(idom,K,L,NY,NX)      = AZERO1(micstt%DOM(idom,K),1.e-11_r8)
      DOM_MicP_drib_vr(idom,K,L,NY,NX) = micstt%DOM_MicP_drib(idom,K)
      SorbedOM_vr(idom,K,L,NY,NX) = AZERO1(micstt%SorbedOM(idom,K),1.e-11_r8)

    enddo
  enddo
  OMBioResdu_vr(1:NumPlantChemElms,1:ndbiomcp,1:KL,L,NY,NX)           = micstt%OMBioResdu(1:NumPlantChemElms,1:ndbiomcp,1:KL)
  mBiomeHeter_vr(1:NumPlantChemElms,1:NumLiveHeterBioms,1:KL,L,NY,NX) = micstt%mBiomeHeter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:KL)
  mBiomeAutor_vr(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)       = micstt%mBiomeAutor(1:NumPlantChemElms,1:NumLiveAutoBioms)
  tRDIM2DOM_col(1:NumPlantChemElms,NY,NX)                             = tRDIM2DOM_col(1:NumPlantChemElms,NY,NX)+micflx%TRDOE2DIE(1:NumPlantChemElms)

  end subroutine MicAPIRecv
end module MicBGCAPI
