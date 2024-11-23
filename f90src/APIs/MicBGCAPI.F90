module MicBGCAPI

  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use SoilBGCNLayMod,            only: DownwardMixOM, sumMicBiomLayL
  use SoilDisturbMod,      only: SOMRemovalByDisturbance
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicrobeDiagTypes,     only: Cumlate_Flux_Diag_type
  use MicForcTypeMod,       only: micforctype
  use minimathmod,          only: AZMAX1
  use EcoSiMParDataMod,     only: micpar
  use MicBGCMod,            only: SoilBGCOneLayer
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
      D998: DO L=0,NL(NY,NX)
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN

             call sumMicBiomLayL(L,NY,NX,OrGM_beg)
             call MicBGC1Layer(I,J,L,NY,NX)
             call sumMicBiomLayL(L,NY,NX,dOrGM)
             dOrGM  = dOrGM-OrGM_beg
             tdOrGM = tdOrGM+dOrGM
          ELSE
            trcg_RMicbTransf_vr(idg_beg:idg_NH3-1,L,NY,NX)   = 0.0_r8
            RNutMicbTransf_vr(ids_NH4B:ids_nuts_end,L,NY,NX) = 0.0_r8
            Micb_N2Fixation_vr(L,NY,NX)                      = 0.0_r8
          ENDIF
  !     MIX LITTER C BETWEEN ADJACENT SOIL LAYERS L AND LL

          call DownwardMixOM(I,J,L,NY,NX)

        ELSE
          trcg_RMicbTransf_vr(idg_beg:idg_NH3-1,L,NY,NX)   = 0.0_r8
          RNutMicbTransf_vr(ids_NH4B:ids_nuts_end,L,NY,NX) = 0.0_r8
          Micb_N2Fixation_vr(L,NY,NX)                      = 0.0_r8
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
  type(Cumlate_Flux_Diag_type) :: naqfdiag
  
  call micflx%ZeroOut()

  call MicAPISend(I,J,L,NY,NX,micfor,micstt,micflx)
  
  call SoilBGCOneLayer(I,J,micfor,micstt,micflx,naqfdiag)

  call MicAPIRecv(I,J,L,NY,NX,micfor%litrm,micstt,micflx,naqfdiag)

  end subroutine Micbgc1Layer
!------------------------------------------------------------------------------------------

  subroutine MicAPISend(I,J,L,NY,NX,micfor,micstt,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: L,NY,NX
  type(micforctype), intent(inout) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx

  integer :: NumMicbFunGrupsPerCmplx, jcplx, k_POM, k_humus
  integer :: ndbiomcp, nlbiomcp, NumMicrobAutrophCmplx, NumHetetrMicCmplx,NE

  NumMicbFunGrupsPerCmplx = micpar%NumMicbFunGrupsPerCmplx
  jcplx                   = micpar%jcplx

  ndbiomcp              = micpar%ndbiomcp
  nlbiomcp              = micpar%nlbiomcp
  NumMicrobAutrophCmplx = micpar%NumMicrobAutrophCmplx
  NumHetetrMicCmplx     = micpar%NumHetetrMicCmplx
  k_humus               = micpar%k_humus
  k_POM                 = micpar%k_POM

  micfor%VOLW0               = VLWatMicP_vr(0,NY,NX)
  micfor%ZERO                = ZERO
  micfor%CCH4E               = AtmGasCgperm3(idg_CH4,NY,NX)
  micfor%COXYE               = AtmGasCgperm3(idg_O2,NY,NX)
  micfor%O2_irrig_conc       = trcVolatile_irrig_conc(idg_O2,NY,NX)
  micfor%O2_rain_conc        = trcVolatile_rain_conc(idg_O2,NY,NX)
  micfor%Irrig2LitRSurf      = Irrig2LitRSurf(NY,NX)
  micfor%Rain2LitRSurf       = Rain2LitRSurf_col(NY,NX)
  micfor%TempOffset          = TempOffset_col(NY,NX)
  micfor%VLitR               = VLitR_col(NY,NX)
  micfor%VWatLitRHoldCapcity = VWatLitRHoldCapcity_col(NY,NX)
  micfor%ZEROS2              = ZEROS2(NY,NX)
  micfor%ZEROS               = ZEROS(NY,NX)
  micfor%VLSoilMicP          = VLSoilMicP_vr(L,NY,NX)
  micfor%THETY               = THETY_vr(L,NY,NX)
  micfor%POROS               = POROS_vr(L,NY,NX)
  micfor%FieldCapacity       = FieldCapacity_vr(L,NY,NX)
  micfor%TKS                 = TKS_vr(L,NY,NX)
  micfor%THETW               = THETW_vr(L,NY,NX)
  micfor%PH                  = PH(L,NY,NX)
  micfor%SoilMicPMassLayer            = VLSoilMicPMass_vr(L,NY,NX)
  micfor%VLSoilPoreMicP               = VLSoilPoreMicP_vr(L,NY,NX)
  micfor%TScal4Difsvity               = TScal4Difsvity_vr(L,NY,NX)
  micfor%VLNOB                        = trcs_VLN_vr(ids_NO3B,L,NY,NX)
  micfor%VLNO3                        = trcs_VLN_vr(ids_NO3,L,NY,NX)
  micfor%VLNH4                        = trcs_VLN_vr(ids_NH4,L,NY,NX)
  micfor%VLNHB                        = trcs_VLN_vr(ids_NH4B,L,NY,NX)
  micfor%VLPO4                        = trcs_VLN_vr(ids_H1PO4,L,NY,NX)
  micfor%VLPOB                        = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
  micfor%PSISoilMatricP               = PSISoilMatricP_vr(L,NY,NX)
  micfor%O2AquaDiffusvity             = SoluteDifusvty_vr(idg_O2,L,NY,NX)
  micfor%ORGC                         = SoilOrgM_vr(ielmc,L,NY,NX)
  micfor%RNO2EcoUptkSoilPrev          = RNO2EcoUptkSoilPrev_vr(L,NY,NX)
  micfor%RN2OEcoUptkSoilPrev          = RN2OEcoUptkSoilPrev_vr(L,NY,NX)
  micfor%RNO2EcoUptkBandPrev          = RNO2EcoUptkBandPrev_vr(L,NY,NX)
  micfor%RO2EcoDmndPrev               = RO2EcoDmndPrev_vr(L,NY,NX)
  micfor%RO2GasXchangePrev            = RO2GasXchangePrev_vr(L,NY,NX)
  micfor%RCH4PhysexchPrev_vr          = RCH4PhysexchPrev_vr(L,NY,NX)
  micfor%RNH4EcoDmndBandPrev          = RNH4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RNO3EcoDmndBandPrev          = RNO3EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RH2PO4EcoDmndBandPrev        = RH2PO4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RH1PO4EcoDmndBandPrev        = RH1PO4EcoDmndBandPrev_vr(L,NY,NX)
  micfor%RO2AquaXchangePrev           = RO2AquaSourcePrev_vr(L,NY,NX)
  micfor%RDOMEcoDmndPrev(1:jcplx)     = RDOMEcoDmndPrev_vr(1:jcplx,L,NY,NX)

!  write(201,*)I+J/24.,L,RDOMEcoDmndPrev_vr(1:jcplx,L,NY,NX)
  micfor%RAcetateEcoDmndPrev(1:jcplx) = RAcetateEcoDmndPrev_vr(1:jcplx,L,NY,NX)
  !is it litter layer?
  micfor%litrm=(L==0)
  !is it surface layer
  micfor%Lsurf=(L==NU(NY,NX))

  !litter layer is modeled
  if(micfor%litrm)then
    micstt%ZNH4TU = AZMAX1(trc_solml_vr(ids_NH4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_NH4B,NU(NY,NX),NY,NX))
    micstt%ZNO3TU = AZMAX1(trc_solml_vr(ids_NO3,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_NO3B,NU(NY,NX),NY,NX))
    micstt%H1P4TU = AZMAX1(trc_solml_vr(ids_H1PO4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_H1PO4B,NU(NY,NX),NY,NX))
    micstt%H2P4TU = AZMAX1(trc_solml_vr(ids_H2PO4,NU(NY,NX),NY,NX))+AZMAX1(trc_solml_vr(ids_H2PO4B,NU(NY,NX),NY,NX))

    micstt%CNH4BU  = trc_solcl_vr(ids_NH4B,NU(NY,NX),NY,NX)
    micstt%CNH4SU  = trc_solcl_vr(ids_NH4,NU(NY,NX),NY,NX)
    micstt%CH2P4U  = trc_solcl_vr(ids_H2PO4,NU(NY,NX),NY,NX)
    micstt%CH2P4BU = trc_solcl_vr(ids_H2PO4B,NU(NY,NX),NY,NX)
    micstt%CH1P4U  = trc_solcl_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    micstt%CH1P4BU = trc_solcl_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
    micstt%CNO3SU  = trc_solcl_vr(ids_NO3,NU(NY,NX),NY,NX)
    micstt%CNO3BU  = trc_solcl_vr(ids_NO3B,NU(NY,NX),NY,NX)

    DO NE=1,NumPlantChemElms
      micstt%SOMPomProtein(NE)  = SolidOM_vr(NE,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)
      micstt%SOMHumProtein(NE)  = SolidOM_vr(NE,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)
      micstt%SOMHumCarbohyd(NE) = SolidOM_vr(NE,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX)
    ENDDO
    micfor%RNH4EcoDmndLitrPrev       = RNH4EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%RNO3EcoDmndLitrPrev       = RNO3EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%RH2PO4EcoDmndLitrPrev     = RH2PO4EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%RH1PO4EcoDmndLitrPrev     = RH1PO4EcoDmndSoilPrev_vr(NU(NY,NX),NY,NX)
    micfor%VOLWU                     = VLWatMicP_vr(NU(NY,NX),NY,NX)
    micfor%ElmAllocmatMicrblitr2POMU = ElmAllocmatMicrblitr2POM_vr(1:2,NU(NY,NX),NY,NX)
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
  micfor%DiffusivitySolutEff(1:NPH) = DiffusivitySolutEff(1:NPH,L,NY,NX)
  micfor%FILM(1:NPH)                = FILM(1:NPH,L,NY,NX)
  micfor%THETPM(1:NPH)              = THETPM(1:NPH,L,NY,NX)
  micfor%VLWatMicPM(1:NPH)          = VLWatMicPM_vr(1:NPH,L,NY,NX)
  micfor%TortMicPM(1:NPH)           = TortMicPM_vr(1:NPH,L,NY,NX)
  micfor%VLsoiAirPM(1:NPH)          = VLsoiAirPM(1:NPH,L,NY,NX)
  micfor%VLsoiAirP                  = VLsoiAirP_vr(L,NY,NX)
  micstt%EPOC                       = EPOC(L,NY,NX)
  micstt%EHUM                       = EHUM(L,NY,NX)
  micstt%ZNH4B             = AZMAX1(trc_solml_vr(ids_NH4B,L,NY,NX))
  micstt%ZNH4S             = AZMAX1(trc_solml_vr(ids_NH4,L,NY,NX))
  micstt%ZNO3B             = AZMAX1(trc_solml_vr(ids_NO3B,L,NY,NX))
  micstt%ZNO3S             = AZMAX1(trc_solml_vr(ids_NO3,L,NY,NX))
  micstt%H1POB             = AZMAX1(trc_solml_vr(ids_H1PO4B,L,NY,NX))
  micstt%H1PO4             = AZMAX1(trc_solml_vr(ids_H1PO4,L,NY,NX))
  micstt%ZNO2B             = AZMAX1(trc_solml_vr(ids_NO2B,L,NY,NX))
  micstt%ZNO2S             = AZMAX1(trc_solml_vr(ids_NO2,L,NY,NX))
  micstt%H2POB             = AZMAX1(trc_solml_vr(ids_H2PO4B,L,NY,NX))
  micstt%H2PO4             = AZMAX1(trc_solml_vr(ids_H2PO4,L,NY,NX))
  micstt%Z2OS              = AZMAX1(trc_solml_vr(idg_N2O,L,NY,NX))  
  micstt%CH4S              = AZMAX1(trc_solml_vr(idg_CH4,L,NY,NX))  
  micstt%OXYS              = AZMAX1(trc_solml_vr(idg_O2,L,NY,NX))  
  micstt%H2GS              = AZMAX1(trc_solml_vr(idg_H2,L,NY,NX))
  micstt%CCO2S             = AZMAX1(trc_solcl_vr(idg_CO2,L,NY,NX))
  micstt%CNO2S             = AZMAX1(trc_solcl_vr(ids_NO2,L,NY,NX))
  micstt%CNO2B             = AZMAX1(trc_solcl_vr(ids_NO2B,L,NY,NX))
  micstt%CZ2OS             = AZMAX1(trc_solcl_vr(idg_N2O,L,NY,NX))
  micstt%COXYS             = AZMAX1(trc_solcl_vr(idg_O2,L,NY,NX))
  micstt%COXYG             = AZMAX1(trc_gascl_vr(idg_O2,L,NY,NX))
  micstt%CZ2GS             = AZMAX1(trc_solcl_vr(idg_N2,L,NY,NX))
  micstt%CH2GS             = AZMAX1(trc_solcl_vr(idg_H2,L,NY,NX))
  micstt%CCH4G             = AZMAX1(trc_gascl_vr(idg_CH4,L,NY,NX))
!  write(115,*)I+J/24.,L,micstt%COXYG,micstt%COXYS,VLsoiAirPM(1,L,NY,NX)
  micstt%O2GSolubility     = GasSolbility_vr(idg_O2,L,NY,NX)  
  micstt%CH4AquaSolubility = GasSolbility_vr(idg_CH4,L,NY,NX)
  micstt%ZNFN0             = ZNFN0(L,NY,NX)
  micstt%ZNFNI             = ZNFNI(L,NY,NX)
  
  micstt%DOM(idom_beg:idom_end,1:jcplx)                    = DOM_vr(idom_beg:idom_end,1:jcplx,L,NY,NX)
  micstt%SorbedOM(idom_beg:idom_end,1:jcplx)               = SorbedOM_vr(idom_beg:idom_end,1:jcplx,L,NY,NX)
  micstt%SolidOMAct(1:jsken,1:jcplx)                       = SolidOMAct_vr(1:jsken,1:jcplx,L,NY,NX)
  micstt%SolidOM(1:NumPlantChemElms,1:jsken,1:jcplx)       = SolidOM_vr(1:NumPlantChemElms,1:jsken,1:jcplx,L,NY,NX)
  micstt%OMBioResdu(1:NumPlantChemElms,1:ndbiomcp,1:jcplx) = OMBioResdu_vr(1:NumPlantChemElms,1:ndbiomcp,1:jcplx,L,NY,NX)
  micstt%CNOSC(1:jsken,1:jcplx)                            = CNOSC(1:jsken,1:jcplx,L,NY,NX)
  micstt%CPOSC(1:jsken,1:jcplx)                            = CPOSC(1:jsken,1:jcplx,L,NY,NX)
  micstt%mBiomeHeter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx)=mBiomeHeter_vr(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,L,NY,NX)
  micstt%mBiomeAutor(1:NumPlantChemElms,1:NumLiveAutoBioms)=mBiomeAutor_vr(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)
  if(.not.micfor%litrm)then
    micfor%AEC  = AEC_vr(L,NY,NX)
    micstt%OXYG = trc_gasml_vr(idg_O2,L,NY,NX)
  endif
  micflx%RNO2DmndSoilChemo=RNO2DmndSoilChemo_vr(L,NY,NX)
  micflx%RNO2DmndBandChemo=RNO2DmndBandChemo_vr(L,NY,NX)
  micflx%RNH4UptkSoilAutor(1:NumMicrobAutrophCmplx)   = RNH4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNH4UptkBandAutor(1:NumMicrobAutrophCmplx)   = RNH4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO3UptkSoilAutor(1:NumMicrobAutrophCmplx)   = RNO3UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RNO3UptkBandAutor(1:NumMicrobAutrophCmplx)   = RNO3UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH2PO4UptkSoilAutor(1:NumMicrobAutrophCmplx) = RH2PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH2PO4UptkBandAutor(1:NumMicrobAutrophCmplx) = RH2PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH1PO4UptkSoilAutor(1:NumMicrobAutrophCmplx) = RH1PO4UptkSoilAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RH1PO4UptkBandAutor(1:NumMicrobAutrophCmplx) = RH1PO4UptkBandAutor_vr(1:NumMicrobAutrophCmplx,L,NY,NX)
  micflx%RO2DmndAutort(1:NumMicrobAutrophCmplx)       = RO2DmndAutort_vr(1:NumMicrobAutrophCmplx,L,NY,NX)

  micflx%RNH4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)   = RNH4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RNH4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)   = RNH4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RNO3DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)   = RNO3DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RNO3DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)   = RNO3DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH2PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx) = RH2PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH2PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx) = RH2PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH1PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx) = RH1PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RH1PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx) = RH1PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx)       = RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RDOCUptkHeter(1:NumHetetrMicCmplx,1:jcplx)       = RDOCUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)
  micflx%RAcetateUptkHeter(1:NumHetetrMicCmplx,1:jcplx)   = RAcetateUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)

  end subroutine MicAPISend

!------------------------------------------------------------------------------------------


  subroutine MicAPIRecv(I,J,L,NY,NX,litrM,micstt,micflx,naqfdiag)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: L,NY,NX
  logical, intent(in) :: litrM
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(in) :: micflx
  type(Cumlate_Flux_Diag_type), intent(in) :: naqfdiag

  integer :: NumMicbFunGrupsPerCmplx, jcplx, NumMicrobAutrophCmplx
  integer :: NE,idom,K
  
  NumMicrobAutrophCmplx = micpar%NumMicrobAutrophCmplx
  NumMicbFunGrupsPerCmplx=micpar%NumMicbFunGrupsPerCmplx
  jcplx=micpar%jcplx

  RCH4ProdHydrog_vr(L,NY,NX) = naqfdiag%tCH4ProdAceto
  RCH4ProdAcetcl_vr(L,NY,NX) = naqfdiag%tCH4ProdH2
  RCH4Oxi_aero_vr(L,NY,NX)   = naqfdiag%tCH4OxiAero
  RFerment_vr(L,NY,NX)       = naqfdiag%tCResp4H2Prod
  RNH3oxi_vr(L,NY,NX)        = naqfdiag%tRNH3Oxi
  RN2ODeniProd_vr(L,NY,NX)   = naqfdiag%TDeniReduxNO2Soil+naqfdiag%TDeniReduxNO2Band
  RN2OChemoProd_vr(L,NY,NX)  = naqfdiag%RN2OProdSoilChemo+naqfdiag%RN2OProdBandChemo
  RN2ONitProd_vr(L,NY,NX)    = naqfdiag%TNitReduxNO2Soil+naqfdiag%TNitReduxNO2Band
  RN2ORedux_vr(L,NY,NX)      = naqfdiag%TReduxN2O

  trcg_RMicbTransf_vr(idg_CO2,L,NY,NX)  = micflx%RCO2NetUptkMicb
  trcg_RMicbTransf_vr(idg_CH4,L,NY,NX)  = micflx%RCH4UptkAutor
  trcg_RMicbTransf_vr(idg_H2,L,NY,NX)   = micflx%RH2NetUptkMicb
  trcg_RMicbTransf_vr(idg_O2,L,NY,NX)   = micflx%RO2UptkMicb
  trcg_RMicbTransf_vr(idg_N2,L,NY,NX)   = micflx%RN2NetUptkMicb
  trcg_RMicbTransf_vr(idg_N2O,L,NY,NX)  = micflx%RN2ONetUptkMicb
  RNutMicbTransf_vr(ids_NH4,L,NY,NX)    = micflx%RNH4MicbTransfSoil
  RNutMicbTransf_vr(ids_NO3,L,NY,NX)    = micflx%RNO3MicbTransfSoil
  RNutMicbTransf_vr(ids_NO2,L,NY,NX)    = micflx%RNO2MicbTransfSoil
  RNutMicbTransf_vr(ids_H2PO4,L,NY,NX)  = micflx%RH2PO4MicbTransfSoil
  RNutMicbTransf_vr(ids_H1PO4,L,NY,NX)  = micflx%RH1PO4MicbTransfSoil
  RNutMicbTransf_vr(ids_NH4B,L,NY,NX)   = micflx%RNH4MicbTransfBand
  RNutMicbTransf_vr(ids_NO3B,L,NY,NX)   = micflx%RNO3MicbTransfBand
  RNutMicbTransf_vr(ids_NO2B,L,NY,NX)   = micflx%RNO2MicbTransfBand
  RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX) = micflx%RH2PO4MicbTransfBand
  RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX) = micflx%RH1PO4MicbTransfBand
  Micb_N2Fixation_vr(L,NY,NX)           = micflx%MicrbN2Fix
  RNO2DmndSoilChemo_vr(L,NY,NX)         = micflx%RNO2DmndSoilChemo
  RNO2DmndBandChemo_vr(L,NY,NX)         = micflx%RNO2DmndBandChemo
  NetNH4Mineralize_CumYr_col(NY,NX)     = NetNH4Mineralize_CumYr_col(NY,NX)+micflx%NetNH4Mineralize
  NetPO4Mineralize_CumYr_col(NY,NX)     = NetPO4Mineralize_CumYr_col(NY,NX)+micflx%NetPO4Mineralize
  DO idom=idom_beg,idom_end
    REcoDOMProd_vr(idom,1:jcplx,L,NY,NX)=micflx%REcoDOMProd(idom,1:jcplx)    
  ENDDO
  RDOMMicProd_vr(idom_beg:idom_end,1:jcplx,L,NY,NX)=REcoDOMProd_vr(idom_beg:idom_end,1:jcplx,L,NY,NX)
 
  RO2DmndAutort_vr(1:NumMicrobAutrophCmplx,L,NY,NX)              = micflx%RO2DmndAutort(1:NumMicrobAutrophCmplx)
  RNH3OxidAutor(1:NumMicrobAutrophCmplx,L,NY,NX)                 = micflx%RNH3OxidAutor(1:NumMicrobAutrophCmplx)
  RNH3OxidAutorBand(1:NumMicrobAutrophCmplx,L,NY,NX)             = micflx%RNH3OxidAutorBand(1:NumMicrobAutrophCmplx)
  RNO2OxidAutor(1:NumMicrobAutrophCmplx,L,NY,NX)                 = micflx%RNO2OxidAutor(1:NumMicrobAutrophCmplx)
  RNO2OxidAutorBand(1:NumMicrobAutrophCmplx,L,NY,NX)             = micflx%RNO2OxidAutorBand(1:NumMicrobAutrophCmplx)
  RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)             = micflx%RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx)
  RDOCUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)          = micflx%RDOCUptkHeter(1:NumHetetrMicCmplx,1:jcplx)
  RAcetateUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)      = micflx%RAcetateUptkHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3ReduxDmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX) = micflx%RNO3ReduxDmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3ReduxDmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX) = micflx%RNO3ReduxDmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO2DmndReduxSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX) = micflx%RNO2DmndReduxSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO2DmndReduxBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX) = micflx%RNO2DmndReduxBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RN2ODmndReduxHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)     = micflx%RN2ODmndReduxHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNH4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)      = micflx%RNH4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNH4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)      = micflx%RNH4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)      = micflx%RNO3DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RNO3DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)      = micflx%RNO3DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH2PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)    = micflx%RH2PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH2PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)    = micflx%RH2PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH1PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)    = micflx%RH1PO4DmndSoilHeter(1:NumHetetrMicCmplx,1:jcplx)
  RH1PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)    = micflx%RH1PO4DmndBandHeter(1:NumHetetrMicCmplx,1:jcplx)
  
  SolidOM_vr(1:NumPlantChemElms,1:jsken,1:jcplx,L,NY,NX) = micstt%SolidOM(1:NumPlantChemElms,1:jsken,1:jcplx)
  SolidOMAct_vr(1:jsken,1:jcplx,L,NY,NX)                 = micstt%SolidOMAct(1:jsken,1:jcplx)
  
  if(litrm)then
    RNH4DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX)   = micflx%RNH4DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RNO3DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX)   = micflx%RNO3DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RH2PO4DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX) = micflx%RH2PO4DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RH1PO4DmndLitrHeter_col(1:NumHetetrMicCmplx,1:jcplx,NY,NX) = micflx%RH1PO4DmndLitrHeter(1:NumHetetrMicCmplx,1:jcplx)
    RNH4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)       = micflx%RNH4UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RNO3UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)       = micflx%RNO3UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RH2PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)     = micflx%RH2PO4UptkLitrAutor(1:NumMicrobAutrophCmplx)
    RH1PO4UptkLitrAutor_col(1:NumMicrobAutrophCmplx,NY,NX)     = micflx%RH1PO4UptkLitrAutor(1:NumMicrobAutrophCmplx)

    DO NE=1,NumPlantChemElms
      SolidOM_vr(NE,micpar%iprotein,micpar%k_POM,NU(NY,NX),NY,NX)    = micstt%SOMPomProtein(NE)
      SolidOM_vr(NE,micpar%iprotein,micpar%k_humus,NU(NY,NX),NY,NX)  = micstt%SOMHumProtein(NE)
      SolidOM_vr(NE,micpar%icarbhyro,micpar%k_humus,NU(NY,NX),NY,NX) = micstt%SOMHumCarbohyd(NE)
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
  RO2UptkSoilM_vr(1:NPH,L,NY,NX)                          = micflx%RO2UptkSoilM(1:NPH)

  TSens4MicbGrwoth_vr(L,NY,NX)                   = micstt%TSens4MicbGrwoth
  VWatMicrobAct_vr(L,NY,NX)                      = micstt%VWatMicrobAct
  TMicHeterAct_vr(L,NY,NX)                       = micstt%TMicHeterAct
  ZNFNI(L,NY,NX)                                 = micstt%ZNFNI
  FracBulkSOMC_vr(1:jcplx,L,NY,NX)               = micstt%FracBulkSOMC(1:jcplx)
  DOM_vr(idom_beg:idom_end,1:jcplx,L,NY,NX)      = micstt%DOM(idom_beg:idom_end,1:jcplx)
  SorbedOM_vr(idom_beg:idom_end,1:jcplx,L,NY,NX) = micstt%SorbedOM(idom_beg:idom_end,1:jcplx)

  do k=1,jcplx
    DO idom=idom_beg,idom_end
      if(abs(SorbedOM_vr(idom,K,L,NY,NX))<1.E-12_r8)SorbedOM_vr(idom,K,L,NY,NX)=0._r8
      if(abs(DOM_vr(idom,K,L,NY,NX))<1.e-12_r8)DOM_vr(idom,K,L,NY,NX)=0._r8
    enddo
  enddo
  OMBioResdu_vr(1:NumPlantChemElms,1:ndbiomcp,1:jcplx,L,NY,NX)           = micstt%OMBioResdu(1:NumPlantChemElms,1:ndbiomcp,1:jcplx)
  mBiomeHeter_vr(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,L,NY,NX) = micstt%mBiomeHeter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx)
  mBiomeAutor_vr(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)          = micstt%mBiomeAutor(1:NumPlantChemElms,1:NumLiveAutoBioms)
  tRDIM2DOM_col(1:NumPlantChemElms,NY,NX)                                = tRDIM2DOM_col(1:NumPlantChemElms,NY,NX)+micflx%TRDOE2DIE(1:NumPlantChemElms)

  end subroutine MicAPIRecv
end module MicBGCAPI
