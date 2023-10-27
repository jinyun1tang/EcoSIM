module PlantAPI
!
! interface to integrate the plant model
  use data_kind_mod   , only : r8 => DAT_KIND_R8
  use ExtractsMod     , only : extracts
  use grosubsMod      , only : grosubs
  use HfuncsMod       , only : hfuncs
  use UptakesMod      , only : uptakes
  use EcoSiMParDataMod, only : micpar, pltpar
  use PlantDisturbMod , only : PrepLandscapeGrazing
  use timings         , only : start_timer, end_timer
  use SoilPhysDataType, only : ALBX
  use EcoSIMSolverPar
  use EcoSIMHistMod
  use SnowDataType
  use TracerIDMod
  use SurfLitterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use EcoSimSumDataType
  use SoilHeatDataType
  use SOMDataType
  use ClimForcDataType
  use EcoSIMCtrlDataType
  use GridDataType
  use RootDataType
  use SoilWaterDataType
  use CanopyDataType
  use PlantDataRateType
  use PlantTraitDataType
  use CanopyRadDataType
  use FlagDataType
  use EcosimBGCFluxType
  use FertilizerDataType
  use SoilBGCDataType
  use PlantMngmtDataType
  use PlantAPIData
implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantModel
  public :: PlantCanopyRadsModel

  contains


  subroutine PlantModel(I,J,NHW,NHE,NVN,NVS)


  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: t1
  integer :: NY,NX

333   FORMAT(A8)

  call PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)

  plt_site%PlantElemntStoreLandscape(:)=PlantElemntStoreLandscape(:)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
!
      call  PlantAPISend(I,J,NY,NX)
!   UPDATE PLANT PHENOLOGY IN 'HFUNC'
!
    !if(lverb)WRITE(*,333)'HFUNC'
    !call start_timer(t1)

      CALL HFUNCs(I,J)
    !call end_timer('HFUNC',t1)

      CALL UPTAKES(I,J)

      CALL GROSUBs(I,J)

      CALL EXTRACTs(I,J)

      call PlantAPIRecv(I,J,NY,NX)
    ENDDO
  ENDDO
  PlantElemntStoreLandscape(:)=plt_site%PlantElemntStoreLandscape(:)


  end subroutine PlantModel
!------------------------------------------------------------------------------------------

  subroutine PlantCanopyRadsModel(I,J,NY,NX,DPTH0)
  use CanopyCondsMod
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DPTH0

  call PlantAPICanMSend(NY,NX)

  call CanopyConditionModel(I,J,DPTH0)

  call PlantAPICanMRecv(NY,NX)

  end subroutine PlantCanopyRadsModel
!------------------------------------------------------------------------------------------

  subroutine PlantAPIRecv(I,J,NY,NX)
  !
  !DESCRIPTION
  !
  use EcoSIMConfig, only : jsken=>jskenc,jcplx => jcplxc
  use PlantAPIData, only : plt_rad
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: NB,NR,NZ,K,L,M,N,I1,NE

  I1=I+1;if(I1>LYRC)I1=1
  IFLGT(NY,NX)=plt_site%IFLGT
  PPT(NY,NX) =plt_site%PPT
  RECO(NY,NX)=plt_bgcr%RECO
  Eco_NBP_col(NY,NX)=plt_bgcr%Eco_NBP_col
  Eco_GPP_col(NY,NX)=plt_bgcr%Eco_GPP_col
  TLEC(NY,NX)=plt_ew%TLEC
  TSHC(NY,NX)=plt_ew%TSHC
  Eco_AutoR_col(NY,NX)=plt_bgcr%Eco_AutoR_col
  ZESNC(1:NumOfPlantChemElements,NY,NX) =plt_bgcr%ZESNC(1:NumOfPlantChemElements)
  XHVSTE(1:NumOfPlantChemElements,NY,NX)=plt_distb%XHVSTE(1:NumOfPlantChemElements)
  CanH2OHeldVg(NY,NX)=plt_ew%CanH2OHeldVg
  TSH(NY,NX)   =plt_ew%TSH
  WTSTGET(1:NumOfPlantChemElements,NY,NX)=plt_biom%WTSTGET(1:NumOfPlantChemElements)
  UVOLO(NY,NX) =plt_ew%UVOLO
  StemAreag(NY,NX) =plt_morph%StemAreag
  CanopyLA_grd(NY,NX) =plt_morph%CanopyLA_grd
  TRN(NY,NX)   =plt_rad%TRN
  TLE(NY,NX)   =plt_ew%TLE
  TGH(NY,NX)   =plt_ew%TGH
  TEVAPP(NY,NX)=plt_ew%TEVAPP
  LWRadCanG(NY,NX) =plt_ew%LWRadCanG
  VapXAir2CanG(NY,NX)=plt_ew%VapXAir2CanG
  THFLXC(NY,NX)=plt_ew%THFLXC
  CanWatg(NY,NX)=plt_ew%CanWatg
  TENGYC(NY,NX)=plt_ew%TENGYC
  TRootGasLoss_disturb(idg_beg:idg_end-1,NY,NX) =plt_rbgc%TRootGasLoss_disturb(idg_beg:idg_end-1)
  TCCAN(NY,NX)=plt_bgcr%TCCAN

  FERT(17:19,I1,NY,NX)=plt_distb%FERT(17:19)
  FERT(3,I1,NY,NX) =plt_distb%FERT(3)
  IYTYP(2,I1,NY,NX)=plt_distb%IYTYP
  FWOODE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs) =plt_allom%FWOODE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)
  FWODRE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs) =plt_allom%FWODRE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)
  FWODBE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs) =plt_allom%FWODBE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)
  FWODLE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)=plt_allom%FWODLE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)
  VOLWOU   =plt_site%VOLWOU
  DO L=1,JC
    WGLFT(L,NY,NX)=plt_biom%WGLFT(L)
    CanopyStemA_lyr(L,NY,NX)=plt_morph%CanopyStemA_lyr(L)
    CanopyLAgrid_lyr(L,NY,NX)=plt_morph%CanopyLAgrid_lyr(L)
  ENDDO
  DO L=NU(NY,NX),NL(NY,NX)
    DO K=1,jcplx
      RDOM_micb_flx(idom_doc,K,L,NY,NX)=plt_bgcr%RDOM_micb_flx(idom_doc,K,L)
      RDOM_micb_flx(idom_don,K,L,NY,NX)=plt_bgcr%RDOM_micb_flx(idom_don,K,L)
      RDOM_micb_flx(idom_dop,K,L,NY,NX)=plt_bgcr%RDOM_micb_flx(idom_dop,K,L)
    ENDDO
    DO M=1,NPH
      ROXSK(M,L,NY,NX)=plt_rbgc%ROXSK(M,L)
    ENDDO
  ENDDO

  DO L=0,JZ
    RPOBX(L,NY,NX) =plt_bgcr%RPOBX(L)
    RP1BX(L,NY,NX) =plt_bgcr%RP1BX(L)
    RN3BX(L,NY,NX) =plt_bgcr%RN3BX(L)
    RNHBX(L,NY,NX) =plt_bgcr%RNHBX(L)
    RP14X(L,NY,NX) =plt_bgcr%RP14X(L)
    RPO4X(L,NY,NX) =plt_bgcr%RPO4X(L)
    RNO3X(L,NY,NX) =plt_bgcr%RNO3X(L)
    RNH4X(L,NY,NX) =plt_bgcr%RNH4X(L)
    ROXYX(L,NY,NX) =plt_bgcr%ROXYX(L)
    THeatRootUptake(L,NY,NX) =plt_ew%THeatRootUptake(L)
    GridPlantRootH2OUptake_vr(L,NY,NX)=plt_ew%GridPlantRootH2OUptake_vr(L)
    DO  K=1,micpar%NumOfPlantLitrCmplxs
      DO NE=1,NumOfPlantChemElements
        DO  M=1,jsken
          LitrfalChemElemnts_vr(NE,M,K,L,NY,NX)=plt_bgcr%LitrfalChemElemnts_vr(NE,M,K,L)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,JZ
    TUPNF(L,NY,NX) =plt_rbgc%TUPNF(L)
    RTDNT(L,NY,NX) =plt_morph%RTDNT(L)
    trcg_TLP(idg_beg:idg_end-1,L,NY,NX)=plt_rbgc%trcg_TLP(idg_beg:idg_end-1,L)
    trcg_TFLA(idg_beg:idg_end-1,L,NY,NX)=plt_rbgc%trcg_TFLA(idg_beg:idg_end-1,L)
    TCO2P(L,NY,NX) =plt_bgcr%TCO2P(L)
    TUPOXP(L,NY,NX)=plt_bgcr%TUPOXP(L)
    TCO2S(L,NY,NX) =plt_bgcr%TCO2S(L)
    TUPOXS(L,NY,NX)=plt_bgcr%TUPOXS(L)
    TUPCHS(L,NY,NX)=plt_bgcr%TUPCHS(L)
    TUPN2S(L,NY,NX)=plt_bgcr%TUPN2S(L)
    TUPN3S(L,NY,NX)=plt_bgcr%TUPN3S(L)
    TUPN3B(L,NY,NX)=plt_bgcr%TUPN3B(L)
    TUPHGS(L,NY,NX)=plt_bgcr%TUPHGS(L)
    TUPNH4(L,NY,NX)=plt_bgcr%TUPNH4(L)
    TUPNO3(L,NY,NX)=plt_bgcr%TUPNO3(L)
    TUPH1B(L,NY,NX)=plt_bgcr%TUPH1B(L)
    TUPH2B(L,NY,NX)=plt_bgcr%TUPH2B(L)
    TUPNHB(L,NY,NX)=plt_bgcr%TUPNHB(L)
    TUPNOB(L,NY,NX)=plt_bgcr%TUPNOB(L)
    TUPH1P(L,NY,NX)=plt_bgcr%TUPH1P(L)
    TUPH2P(L,NY,NX)=plt_bgcr%TUPH2P(L)
    DO  K=1,jcplx
      TDFOME(1:NumOfPlantChemElements,K,L,NY,NX)=plt_bgcr%TDFOME(1:NumOfPlantChemElements,K,L)
    ENDDO
  ENDDO
  DO NZ=1,NP0(NY,NX)
    WTRTE(1:NumOfPlantChemElements,NZ,NY,NX) =plt_biom%WTRTE(1:NumOfPlantChemElements,NZ)
    BALE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_site%BALE(1:NumOfPlantChemElements,NZ)
    CanopyNonstructElements_pft(1:NumOfPlantChemElements,NZ,NY,NX)=plt_biom%CanopyNonstructElements_pft(1:NumOfPlantChemElements,NZ)
    EPOLNP(1:NumOfPlantChemElements,NZ,NY,NX)=plt_biom%EPOLNP(1:NumOfPlantChemElements,NZ)
    CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ,NY,NX)=plt_biom%CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ)
    HESNC(1:NumOfPlantChemElements,NZ,NY,NX) =plt_bgcr%HESNC(1:NumOfPlantChemElements,NZ)
    HVSTE(1:NumOfPlantChemElements,NZ,NY,NX) =plt_distb%HVSTE(1:NumOfPlantChemElements,NZ)
    RSETE(1:NumOfPlantChemElements,NZ,NY,NX) =plt_pheno%RSETE(1:NumOfPlantChemElements,NZ)
    TESN0(1:NumOfPlantChemElements,NZ,NY,NX) =plt_bgcr%TESN0(1:NumOfPlantChemElements,NZ)
    TESNC(1:NumOfPlantChemElements,NZ,NY,NX) =plt_bgcr%TESNC(1:NumOfPlantChemElements,NZ)
    THVSTE(1:NumOfPlantChemElements,NZ,NY,NX) =plt_distb%THVSTE(1:NumOfPlantChemElements,NZ)
    TEUPTK(1:NumOfPlantChemElements,NZ,NY,NX) =plt_rbgc%TEUPTK(1:NumOfPlantChemElements,NZ)
    UPOME(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_rbgc%UPOME(1:NumOfPlantChemElements,NZ)
    WTSTGE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTSTGE(1:NumOfPlantChemElements,NZ)
    WTRVE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTRVE(1:NumOfPlantChemElements,NZ)
    CanPShootElmMass(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%CanPShootElmMass(1:NumOfPlantChemElements,NZ)
    WTLFE(1:NumOfPlantChemElements,NZ,NY,NX)   =plt_biom%WTLFE(1:NumOfPlantChemElements,NZ)
    WTSHEE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTSHEE(1:NumOfPlantChemElements,NZ)
    WTSTKE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTSTKE(1:NumOfPlantChemElements,NZ)
    WTRSVE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTRSVE(1:NumOfPlantChemElements,NZ)
    WTHSKE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTHSKE(1:NumOfPlantChemElements,NZ)
    WTEARE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTEARE(1:NumOfPlantChemElements,NZ)
    WTGRE(1:NumOfPlantChemElements,NZ,NY,NX)   =plt_biom%WTGRE(1:NumOfPlantChemElements,NZ)
    WTRTSE(1:NumOfPlantChemElements,NZ,NY,NX)  =plt_biom%WTRTSE(1:NumOfPlantChemElements,NZ)
    WTNDE(1:NumOfPlantChemElements,NZ,NY,NX)   =plt_biom%WTNDE(1:NumOfPlantChemElements,NZ)
    HEUPTK(1:NumOfPlantChemElements,NZ,NY,NX)=plt_rbgc%HEUPTK(1:NumOfPlantChemElements,NZ)
    CanopyLeafA_pft(NZ,NY,NX) =plt_morph%CanopyLeafA_pft(NZ)
    CanPSA(NZ,NY,NX) =plt_morph%CanPSA(NZ)
    NoduleNonstructCconc_pft(NZ,NY,NX)=plt_biom%NoduleNonstructCconc_pft(NZ)
    CO2NetFix_pft(NZ,NY,NX)  =plt_bgcr%CO2NetFix_pft(NZ)
    CO2Q(NZ,NY,NX)  =plt_photo%CO2Q(NZ)
    CO2I(NZ,NY,NX)  =plt_photo%CO2I(NZ)
    CO2L(NZ,NY,NX)  =plt_photo%CO2L(NZ)
    ClumpFactor(NZ,NY,NX)    =plt_morph%ClumpFactor(NZ)
    CNWS(NZ,NY,NX)  =plt_allom%CNWS(NZ)
    CPWS(NZ,NY,NX)  =plt_allom%CPWS(NZ)
    RootFracRemobilizableBiom(NZ,NY,NX) =plt_allom%RootFracRemobilizableBiom(NZ)
    CNRTS(NZ,NY,NX) =plt_allom%CNRTS(NZ)
    CPRTS(NZ,NY,NX) =plt_allom%CPRTS(NZ)
    ETCanP(NZ,NY,NX) =plt_ew%ETCanP(NZ)
    CARBN(NZ,NY,NX) =plt_bgcr%CARBN(NZ)
    CHILL(NZ,NY,NX) =plt_photo%CHILL(NZ)
    DCO2(NZ,NY,NX)  =plt_photo%DCO2(NZ)
    DTKC(NZ,NY,NX)  =plt_ew%DTKC(NZ)
    ENGYX(NZ,NY,NX) =plt_ew%ENGYX(NZ)
    PTrans(NZ,NY,NX)    =plt_ew%PTrans(NZ)
    VapXAir2PCan(NZ,NY,NX) =plt_ew%VapXAir2PCan(NZ)
    EvapTransHeatP(NZ,NY,NX) =plt_ew%EvapTransHeatP(NZ)
    FMOL(NZ,NY,NX)  =plt_photo%FMOL(NZ)
    FracPARByCanP(NZ,NY,NX) =plt_rad%FracPARByCanP(NZ)
    FNOD(NZ,NY,NX)  =plt_allom%FNOD(NZ)
    GRNO(NZ,NY,NX)  =plt_morph%GRNO(NZ)
    HypoctoylHeight(NZ,NY,NX) =plt_morph%HypoctoylHeight(NZ)
    HTC(NZ,NY,NX)   =plt_pheno%HTC(NZ)
    HeatStorCanP(NZ,NY,NX) =plt_ew%HeatStorCanP(NZ)
    CanPHeight4WatUptake(NZ,NY,NX) =plt_morph%CanPHeight4WatUptake(NZ)
    IFLGC(NZ,NY,NX) =plt_pheno%IFLGC(NZ)
    IDTH(NZ,NY,NX)  =plt_pheno%IDTH(NZ)
    IDTHP(NZ,NY,NX) =plt_pheno%IDTHP(NZ)
    IDTHR(NZ,NY,NX) =plt_pheno%IDTHR(NZ)
    IDAY0(NZ,NY,NX) =plt_distb%IDAY0(NZ)
    IYR0(NZ,NY,NX)  =plt_distb%IYR0(NZ)
    IFLGI(NZ,NY,NX) =plt_pheno%IFLGI(NZ)
    IDAYH(NZ,NY,NX) =plt_distb%IDAYH(NZ)
    IYRH(NZ,NY,NX)  =plt_distb%IYRH(NZ)
    NIXBotRootLayer(NZ,NY,NX)   =plt_morph%NIXBotRootLayer(NZ)
    NBT(NZ,NY,NX)   =plt_morph%NBT(NZ)
    NumOfBranches_pft(NZ,NY,NX)   =plt_morph%NumOfBranches_pft(NZ)
    NRT(NZ,NY,NX)   =plt_morph%NRT(NZ)
    NI(NZ,NY,NX)    =plt_morph%NI(NZ)
    NGTopRootLayer(NZ,NY,NX)    =plt_morph%NGTopRootLayer(NZ)
    NB1(NZ,NY,NX)   =plt_morph%NB1(NZ)
    NNOD(NZ,NY,NX)  =plt_morph%NNOD(NZ)
    O2L(NZ,NY,NX)   =plt_photo%O2L(NZ)
    O2I(NZ,NY,NX)   =plt_photo%O2I(NZ)
    OFFST(NZ,NY,NX) =plt_pheno%OFFST(NZ)
    PlantO2Stress(NZ,NY,NX)  =plt_pheno%PlantO2Stress(NZ)
    PPX(NZ,NY,NX)   =plt_site%PPX(NZ)
    pftPlantPopulation(NZ,NY,NX)    =plt_site%pftPlantPopulation(NZ)
    PPI(NZ,NY,NX)   =plt_site%PPI(NZ)
    PSICanP(NZ,NY,NX) =plt_ew%PSICanP(NZ)
    PSICanPOsmo(NZ,NY,NX) =plt_ew%PSICanPOsmo(NZ)
    PSICanPTurg(NZ,NY,NX) =plt_ew%PSICanPTurg(NZ)
    PSICanPDailyMin(NZ,NY,NX) =plt_ew%PSICanPDailyMin(NZ)

    RootGasLoss_disturb(idg_beg:idg_end-1,NZ,NY,NX) =plt_bgcr%RootGasLoss_disturb(idg_beg:idg_end-1,NZ)
    MinCanPStomaResistH2O(NZ,NY,NX)  =plt_photo%MinCanPStomaResistH2O(NZ)
    MaxCanPStomaResistH2O(NZ,NY,NX)  =plt_photo%MaxCanPStomaResistH2O(NZ)
    RCMX(NZ,NY,NX)  =plt_photo%RCMX(NZ)
    RNH3C(NZ,NY,NX) =plt_bgcr%RNH3C(NZ)
    RadNet2CanP(NZ,NY,NX)  =plt_rad%RadNet2CanP(NZ)
    RAZ(NZ,NY,NX)   =plt_ew%RAZ(NZ)
    CanPStomaResistH2O(NZ,NY,NX)    =plt_photo%CanPStomaResistH2O(NZ)
    CanPbndlResist(NZ,NY,NX)    =plt_photo%CanPbndlResist(NZ)
    PlantinDepth(NZ,NY,NX)=plt_morph%PlantinDepth(NZ)
    SCO2(NZ,NY,NX)  =plt_photo%SCO2(NZ)
    SO2(NZ,NY,NX)   =plt_photo%SO2(NZ)
    SSTX(NZ,NY,NX)  =plt_pheno%SSTX(NZ)
    SeedVolume(NZ,NY,NX)  =plt_morph%SeedVolume(NZ)
    SeedLength(NZ,NY,NX)  =plt_morph%SeedLength(NZ)
    SeedArea(NZ,NY,NX)  =plt_morph%SeedArea(NZ)
    SeedinDepth(NZ,NY,NX) =plt_morph%SeedinDepth(NZ)
    HeatXAir2PCan(NZ,NY,NX) =plt_ew%HeatXAir2PCan(NZ)
    TCO2T(NZ,NY,NX) =plt_bgcr%TCO2T(NZ)
    TCO2A(NZ,NY,NX) =plt_bgcr%TCO2A(NZ)
    TCZ(NZ,NY,NX)   =plt_pheno%TCZ(NZ)
    TCX(NZ,NY,NX)   =plt_pheno%TCX(NZ)
    TNH3C(NZ,NY,NX)  =plt_bgcr%TNH3C(NZ)
    TZUPFX(NZ,NY,NX) =plt_bgcr%TZUPFX(NZ)
    TKC(NZ,NY,NX)    =plt_ew%TKC(NZ)
    TCC(NZ,NY,NX)    =plt_ew%TCC(NZ)
    LWRadCanP(NZ,NY,NX)  =plt_rad%LWRadCanP(NZ)
    TKCZ(NZ,NY,NX)   =plt_ew%TKCZ(NZ)
    TKG(NZ,NY,NX)    =plt_pheno%TKG(NZ)
    TCG(NZ,NY,NX)    =plt_pheno%TCG(NZ)
    fTgrowCanP(NZ,NY,NX)   =plt_pheno%fTgrowCanP(NZ)
    UPNF(NZ,NY,NX)   =plt_rbgc%UPNF(NZ)
    UPNH4(NZ,NY,NX)  =plt_rbgc%UPNH4(NZ)
    UPNO3(NZ,NY,NX)  =plt_rbgc%UPNO3(NZ)
    UPH2P(NZ,NY,NX)  =plt_rbgc%UPH2P(NZ)
    UPH1P(NZ,NY,NX)  =plt_rbgc%UPH1P(NZ)
    VHeatCapCanP(NZ,NY,NX)  =plt_ew%VHeatCapCanP(NZ)
    CanWatP(NZ,NY,NX)  =plt_ew%CanWatP(NZ)
    VCO2F(NZ,NY,NX)  =plt_distb%VCO2F(NZ)
    VCH4F(NZ,NY,NX)  =plt_distb%VCH4F(NZ)
    VOXYF(NZ,NY,NX)  =plt_distb%VOXYF(NZ)
    VNH3F(NZ,NY,NX)  =plt_distb%VNH3F(NZ)
    VN2OF(NZ,NY,NX)  =plt_distb%VN2OF(NZ)
    VPO4F(NZ,NY,NX)  =plt_distb%VPO4F(NZ)
    WatByPCan(NZ,NY,NX)  =plt_ew%WatByPCan(NZ)
    WSTR(NZ,NY,NX)   =plt_pheno%WSTR(NZ)
    WTRVX(NZ,NY,NX)  =plt_biom%WTRVX(NZ)
    CanPStalkC(NZ,NY,NX)  =plt_biom%CanPStalkC(NZ)
    CanopyLeafShethC_pft(NZ,NY,NX)   =plt_biom%CanopyLeafShethC_pft(NZ)
    WTRTA(NZ,NY,NX)  =plt_biom%WTRTA(NZ)
    XKCO2L(NZ,NY,NX) =plt_photo%XKCO2L(NZ)
    XKCO2O(NZ,NY,NX) =plt_photo%XKCO2O(NZ)
    CanopyHeight(NZ,NY,NX)     =plt_morph%CanopyHeight(NZ)
    ZNPP(NZ,NY,NX)   =plt_bgcr%ZNPP(NZ)
    ZEROP(NZ,NY,NX)  =plt_biom%ZEROP(NZ)
    ZEROQ(NZ,NY,NX)  =plt_rbgc%ZEROQ(NZ)
    ZEROL(NZ,NY,NX)  =plt_biom%ZEROL(NZ)
    ZTYP(NZ,NY,NX)   =plt_pheno%ZTYP(NZ)
    HVST(NZ,I,NY,NX) =plt_distb%HVST(NZ)
    IHVST(NZ,I,NY,NX)=plt_distb%IHVST(NZ)
    JHVST(NZ,I,NY,NX)=plt_distb%JHVST(NZ)
    THIN_pft(NZ,I,NY,NX) =plt_distb%THIN_pft(NZ)
    DO L=1,JZ
      WTNDLE(1:NumOfPlantChemElements,L,NZ,NY,NX) =plt_biom%WTNDLE(1:NumOfPlantChemElements,L,NZ)
      EPOOLN(1:NumOfPlantChemElements,L,NZ,NY,NX)=plt_biom%EPOOLN(1:NumOfPlantChemElements,L,NZ)
      RUPNF(L,NZ,NY,NX) =plt_bgcr%RUPNF(L,NZ)
      fTgrowRootP(L,NZ,NY,NX)  =plt_pheno%fTgrowRootP(L,NZ)
    ENDDO
    DO L=1,JC
      CanopyLeafApft_lyr(L,NZ,NY,NX)=plt_morph%CanopyLeafApft_lyr(L,NZ)
      CanopyLeafCpft_lyr(L,NZ,NY,NX)=plt_biom%CanopyLeafCpft_lyr(L,NZ)
      CanopyStemApft_lyr(L,NZ,NY,NX)=plt_morph%CanopyStemApft_lyr(L,NZ)
    ENDDO

    DO L=0,JZ
      DO K=1,micpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          ESNC(1:NumOfPlantChemElements,M,K,L,NZ,NY,NX)=plt_bgcr%ESNC(1:NumOfPlantChemElements,M,K,L,NZ)
        ENDDO
      ENDDO
    ENDDO

    
    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO NE=1,NumOfPlantChemElements
        EPOOL(NE,NB,NZ,NY,NX) =plt_biom%EPOOL(NE,NB,NZ)
        WTEARBE(NE,NB,NZ,NY,NX)=plt_biom%WTEARBE(NE,NB,NZ)
      ENDDO
    ENDDO


    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO NE=1,NumOfPlantChemElements
        EPOLNB(NE,NB,NZ,NY,NX)=plt_biom%EPOLNB(NE,NB,NZ)
        WTSHTBE(NE,NB,NZ,NY,NX)=plt_biom%WTSHTBE(NE,NB,NZ)
        WTSHEBE(NE,NB,NZ,NY,NX)=plt_biom%WTSHEBE(NE,NB,NZ)
        WTSTKBE(NE,NB,NZ,NY,NX)=plt_biom%WTSTKBE(NE,NB,NZ)
        CEPOLB(NE,NB,NZ,NY,NX)=plt_biom%CEPOLB(NE,NB,NZ)
        WTLFBE(NE,NB,NZ,NY,NX) =plt_biom%WTLFBE(NE,NB,NZ)
        WTRSVBE(NE,NB,NZ,NY,NX)=plt_biom%WTRSVBE(NE,NB,NZ)
        WTHSKBE(NE,NB,NZ,NY,NX)=plt_biom%WTHSKBE(NE,NB,NZ)
        WTGRBE(NE,NB,NZ,NY,NX) =plt_biom%WTGRBE(NE,NB,NZ)
        WTNDBE(NE,NB,NZ,NY,NX) =plt_biom%WTNDBE(NE,NB,NZ)
        WGSHEXE(NE,NB,NZ,NY,NX)=plt_biom%WGSHEXE(NE,NB,NZ)
      ENDDO
    ENDDO
    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO L=1,JC
        CanopyBranchStemApft_lyr(L,NB,NZ,NY,NX)=plt_morph%CanopyBranchStemApft_lyr(L,NB,NZ)
      ENDDO
      HourCounter4LeafOut_brch(NB,NZ,NY,NX)  =plt_pheno%HourCounter4LeafOut_brch(NB,NZ)
      CanopyBranchLeafA_pft(NB,NZ,NY,NX) =plt_morph%CanopyBranchLeafA_pft(NB,NZ)
      ARLFZ(NB,NZ,NY,NX) =plt_morph%ARLFZ(NB,NZ)

      DGSTGI(NB,NZ,NY,NX)=plt_pheno%DGSTGI(NB,NZ)
      DGSTGF(NB,NZ,NY,NX)=plt_pheno%DGSTGF(NB,NZ)
      FLG4(NB,NZ,NY,NX)  =plt_pheno%FLG4(NB,NZ)
      RubiscoActivity_brpft(NB,NZ,NY,NX)  =plt_photo%RubiscoActivity_brpft(NB,NZ)
      FDBKX(NB,NZ,NY,NX) =plt_photo%FDBKX(NB,NZ)
      FLGZ(NB,NZ,NY,NX)  =plt_pheno%FLGZ(NB,NZ)
      GROUP(NB,NZ,NY,NX) =plt_pheno%GROUP(NB,NZ)
      GSTGI(NB,NZ,NY,NX) =plt_pheno%GSTGI(NB,NZ)
      GSTGF(NB,NZ,NY,NX) =plt_pheno%GSTGF(NB,NZ)
      GRNXB(NB,NZ,NY,NX) =plt_morph%GRNXB(NB,NZ)
      GRNOB(NB,NZ,NY,NX) =plt_morph%GRNOB(NB,NZ)
      GRWTB(NB,NZ,NY,NX) =plt_allom%GRWTB(NB,NZ)
      CanPBranchHeight(NB,NZ,NY,NX)=plt_morph%CanPBranchHeight(NB,NZ)
      IFLGP(NB,NZ,NY,NX) =plt_pheno%IFLGP(NB,NZ)
      IDTHB(NB,NZ,NY,NX) =plt_pheno%IDTHB(NB,NZ)
      IFLGF(NB,NZ,NY,NX) =plt_pheno%IFLGF(NB,NZ)
      IFLGE(NB,NZ,NY,NX) =plt_pheno%IFLGE(NB,NZ)
      IFLGA(NB,NZ,NY,NX) =plt_pheno%IFLGA(NB,NZ)
      IFLGG(NB,NZ,NY,NX) =plt_pheno%IFLGG(NB,NZ)
      IFLGR(NB,NZ,NY,NX) =plt_pheno%IFLGR(NB,NZ)
      IFLGQ(NB,NZ,NY,NX) =plt_pheno%IFLGQ(NB,NZ)
      KVSTG(NB,NZ,NY,NX) =plt_pheno%KVSTG(NB,NZ)
      KLEAF(NB,NZ,NY,NX) =plt_morph%KLEAF(NB,NZ)
      KLEAFX(NB,NZ,NY,NX)=plt_morph%KLEAFX(NB,NZ)
      KVSTGN(NB,NZ,NY,NX)=plt_pheno%KVSTGN(NB,NZ)
      BranchNumber_brchpft(NB,NZ,NY,NX)  =plt_morph%BranchNumber_brchpft(NB,NZ)
      PSTG(NB,NZ,NY,NX)  =plt_morph%PSTG(NB,NZ)
      PSTGI(NB,NZ,NY,NX) =plt_morph%PSTGI(NB,NZ)
      PSTGF(NB,NZ,NY,NX) =plt_morph%PSTGF(NB,NZ)
      RCELX(1:NumOfPlantChemElements,NB,NZ,NY,NX) =plt_pheno%RCELX(1:NumOfPlantChemElements,NB,NZ)
      RCESX(1:NumOfPlantChemElements,NB,NZ,NY,NX) =plt_pheno%RCESX(1:NumOfPlantChemElements,NB,NZ)
      RNH3B(NB,NZ,NY,NX) =plt_rbgc%RNH3B(NB,NZ)
      TGSTGI(NB,NZ,NY,NX)=plt_pheno%TGSTGI(NB,NZ)
      TGSTGF(NB,NZ,NY,NX)=plt_pheno%TGSTGF(NB,NZ)
      VRNY(NB,NZ,NY,NX)  =plt_pheno%VRNY(NB,NZ)
      VRNZ(NB,NZ,NY,NX)  =plt_pheno%VRNZ(NB,NZ)
      VRNS(NB,NZ,NY,NX)  =plt_pheno%VRNS(NB,NZ)
      VRNF(NB,NZ,NY,NX)  =plt_pheno%VRNF(NB,NZ)
      VSTG(NB,NZ,NY,NX)  =plt_morph%VSTG(NB,NZ)
      VSTGX(NB,NZ,NY,NX) =plt_pheno%VSTGX(NB,NZ)
      CanPBLeafShethC(NB,NZ,NY,NX) =plt_biom%CanPBLeafShethC(NB,NZ)

      WGLFEX(1:NumOfPlantChemElements,NB,NZ,NY,NX) =plt_biom%WGLFEX(1:NumOfPlantChemElements,NB,NZ)
      WTSTXBE(1:NumOfPlantChemElements,NB,NZ,NY,NX)=plt_biom%WTSTXBE(1:NumOfPlantChemElements,NB,NZ)
      CanPBStalkC(NB,NZ,NY,NX)=plt_biom%CanPBStalkC(NB,NZ)

      DO K=0,JNODS
        ARLF(K,NB,NZ,NY,NX)=plt_morph%ARLF1(K,NB,NZ)
        HTNODX(K,NB,NZ,NY,NX)=plt_morph%HTNODX(K,NB,NZ)
        HTNODE(K,NB,NZ,NY,NX)=plt_morph%HTNODE(K,NB,NZ)
        CanPSheathHeight(K,NB,NZ,NY,NX) =plt_morph%CanPSheathHeight(K,NB,NZ)
        WGNODE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)=plt_biom%WGNODE(1:NumOfPlantChemElements,K,NB,NZ)
        WGLFE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)  =plt_biom%WGLFE(1:NumOfPlantChemElements,K,NB,NZ)
        WSLF(K,NB,NZ,NY,NX)  =plt_biom%WSLF(K,NB,NZ)
        WGSHE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX) =plt_biom%WGSHE(1:NumOfPlantChemElements,K,NB,NZ)
        WSSHE(K,NB,NZ,NY,NX) =plt_biom%WSSHE(K,NB,NZ)
      ENDDO
      DO  L=1,JC
        DO N=1,JLI
          StemA_lyrnodbrchpft(N,L,NB,NZ,NY,NX)=plt_morph%StemA_lyrnodbrchpft(N,L,NB,NZ)
        ENDDO
      ENDDO
      DO K=0,JNODS
        DO  L=1,JC
          CanPLNBLA(L,K,NB,NZ,NY,NX) =plt_morph%CanPLNBLA(L,K,NB,NZ)
          WGLFLE(1:NumOfPlantChemElements,L,K,NB,NZ,NY,NX) =plt_biom%WGLFLE(1:NumOfPlantChemElements,L,K,NB,NZ)
        ENDDO
      ENDDO
      DO M=1,pltpar%NumGrothStages
        IDAY(M,NB,NZ,NY,NX)=plt_pheno%IDAY(M,NB,NZ)
      ENDDO
      DO K=1,JNODS
        DO  L=1,JC
          DO N=1,JLI
            LeafA_lyrnodbrchpft(N,L,K,NB,NZ,NY,NX) =plt_morph%LeafA_lyrnodbrchpft(N,L,K,NB,NZ)
            LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ,NY,NX)=plt_photo%LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ)
          ENDDO
        ENDDO
        CPOOL3(K,NB,NZ,NY,NX)=plt_photo%CPOOL3(K,NB,NZ)
        CPOOL4(K,NB,NZ,NY,NX)=plt_photo%CPOOL4(K,NB,NZ)
        CO2B(K,NB,NZ,NY,NX)  =plt_photo%CO2B(K,NB,NZ)
        COMPL(K,NB,NZ,NY,NX) =plt_photo%COMPL(K,NB,NZ)
        CBXN(K,NB,NZ,NY,NX)  =plt_photo%CBXN(K,NB,NZ)
        CBXN4(K,NB,NZ,NY,NX) =plt_photo%CBXN4(K,NB,NZ)
        ETGRO(K,NB,NZ,NY,NX) =plt_photo%ETGRO(K,NB,NZ)
        ETGR4(K,NB,NZ,NY,NX) =plt_photo%ETGR4(K,NB,NZ)
        FDBK4(K,NB,NZ,NY,NX) =plt_photo%FDBK4(K,NB,NZ)
        HCOB(K,NB,NZ,NY,NX)  =plt_photo%HCOB(K,NB,NZ)
        VCGRO(K,NB,NZ,NY,NX) =plt_photo%VCGRO(K,NB,NZ)
        VGRO(K,NB,NZ,NY,NX)  =plt_photo%VGRO(K,NB,NZ)
        VCGR4(K,NB,NZ,NY,NX) =plt_photo%VCGR4(K,NB,NZ)
        VGRO4(K,NB,NZ,NY,NX) =plt_photo%VGRO4(K,NB,NZ)
      ENDDO
    ENDDO
    DO M=1,jsken
      WTSTDE(1:NumOfPlantChemElements,M,NZ,NY,NX)=plt_biom%WTSTDE(1:NumOfPlantChemElements,M,NZ)
    ENDDO

    DO  L=1,JZ
      DO N=1,pltpar%jroots
        EPOOLR(1:NumOfPlantChemElements,N,L,NZ,NY,NX)=plt_biom%EPOOLR(1:NumOfPlantChemElements,N,L,NZ)
        RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ,NY,NX)=plt_biom%RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ)
        RootProteinConc_pftvr(N,L,NZ,NY,NX)=plt_biom%RootProteinConc_pftvr(N,L,NZ)
        trcg_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)  =plt_rbgc%trcg_rootml(idg_beg:idg_end-1,N,L,NZ)
        trcs_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)  =plt_rbgc%trcs_rootml(idg_beg:idg_end-1,N,L,NZ)
        PSIRoot(N,L,NZ,NY,NX) =plt_ew%PSIRoot(N,L,NZ)
        PSIRootOSMO(N,L,NZ,NY,NX) =plt_ew%PSIRootOSMO(N,L,NZ)
        PSIRootTurg(N,L,NZ,NY,NX) =plt_ew%PSIRootTurg(N,L,NZ)
        PrimRootXNumL(N,L,NZ,NY,NX)  =plt_morph%PrimRootXNumL(N,L,NZ)
        SecndRootXNumL(N,L,NZ,NY,NX)  =plt_morph%SecndRootXNumL(N,L,NZ)
        RootLenPerP(N,L,NZ,NY,NX) =plt_morph%RootLenPerP(N,L,NZ)
        RootLenDensNLP(N,L,NZ,NY,NX) =plt_morph%RootLenDensNLP(N,L,NZ)
        RTVLP(N,L,NZ,NY,NX) =plt_morph%RTVLP(N,L,NZ)
        RTVLW(N,L,NZ,NY,NX) =plt_morph%RTVLW(N,L,NZ)
        PrimRootRadius(N,L,NZ,NY,NX) =plt_morph%PrimRootRadius(N,L,NZ)
        SecndRootRadius(N,L,NZ,NY,NX) =plt_morph%SecndRootRadius(N,L,NZ)
        RTARP(N,L,NZ,NY,NX) =plt_morph%RTARP(N,L,NZ)
        AveSecndRootLen(N,L,NZ,NY,NX) =plt_morph%AveSecndRootLen(N,L,NZ)
        RCO2M(N,L,NZ,NY,NX) =plt_rbgc%RCO2M(N,L,NZ)
        RCO2N(N,L,NZ,NY,NX) =plt_rbgc%RCO2N(N,L,NZ)
        RCO2A(N,L,NZ,NY,NX) =plt_rbgc%RCO2A(N,L,NZ)
        RCO2P(N,L,NZ,NY,NX) =plt_rbgc%RCO2P(N,L,NZ)
        RUPOXP(N,L,NZ,NY,NX)=plt_rbgc%RUPOXP(N,L,NZ)
        RCO2S(N,L,NZ,NY,NX) =plt_rbgc%RCO2S(N,L,NZ)
        RUPOXS(N,L,NZ,NY,NX)=plt_rbgc%RUPOXS(N,L,NZ)
        RUPCHS(N,L,NZ,NY,NX)=plt_rbgc%RUPCHS(N,L,NZ)
        RUPN2S(N,L,NZ,NY,NX)=plt_rbgc%RUPN2S(N,L,NZ)
        RUPN3S(N,L,NZ,NY,NX)=plt_rbgc%RUPN3S(N,L,NZ)
        RUPN3B(N,L,NZ,NY,NX)=plt_rbgc%RUPN3B(N,L,NZ)
        RUPHGS(N,L,NZ,NY,NX)=plt_rbgc%RUPHGS(N,L,NZ)
        trcg_RFLA(idg_CO2,N,L,NZ,NY,NX)=plt_rbgc%trcg_RFLA(idg_CO2,N,L,NZ)
        trcg_RFLA(idg_O2,N,L,NZ,NY,NX)=plt_rbgc%trcg_RFLA(idg_O2,N,L,NZ)
        trcg_RFLA(idg_CH4,N,L,NZ,NY,NX)=plt_rbgc%trcg_RFLA(idg_CH4,N,L,NZ)
        trcg_RFLA(idg_N2O,N,L,NZ,NY,NX)=plt_rbgc%trcg_RFLA(idg_N2O,N,L,NZ)
        trcg_RFLA(idg_NH3,N,L,NZ,NY,NX)=plt_rbgc%trcg_RFLA(idg_NH3,N,L,NZ)
        trcg_RFLA(idg_H2,N,L,NZ,NY,NX)=plt_rbgc%trcg_RFLA(idg_H2,N,L,NZ)
        trcg_RDFA(idg_CO2,N,L,NZ,NY,NX)=plt_rbgc%trcg_RDFA(idg_CO2,N,L,NZ)
        trcg_RDFA(idg_O2,N,L,NZ,NY,NX)=plt_rbgc%trcg_RDFA(idg_O2,N,L,NZ)
        trcg_RDFA(idg_CH4,N,L,NZ,NY,NX)=plt_rbgc%trcg_RDFA(idg_CH4,N,L,NZ)
        trcg_RDFA(idg_N2O,N,L,NZ,NY,NX)=plt_rbgc%trcg_RDFA(idg_N2O,N,L,NZ)
        trcg_RDFA(idg_NH3,N,L,NZ,NY,NX)=plt_rbgc%trcg_RDFA(idg_NH3,N,L,NZ)
        trcg_RDFA(idg_H2,N,L,NZ,NY,NX)=plt_rbgc%trcg_RDFA(idg_H2,N,L,NZ)
        RUNNHP(N,L,NZ,NY,NX)=plt_rbgc%RUNNHP(N,L,NZ)
        RUPNH4(N,L,NZ,NY,NX)=plt_rbgc%RUPNH4(N,L,NZ)
        RUONH4(N,L,NZ,NY,NX)=plt_rbgc%RUONH4(N,L,NZ)
        RUCNH4(N,L,NZ,NY,NX)=plt_rbgc%RUCNH4(N,L,NZ)
        RUNNBP(N,L,NZ,NY,NX)=plt_rbgc%RUNNBP(N,L,NZ)
        RUPNHB(N,L,NZ,NY,NX)=plt_rbgc%RUPNHB(N,L,NZ)
        RUONHB(N,L,NZ,NY,NX)=plt_rbgc%RUONHB(N,L,NZ)
        RUCNHB(N,L,NZ,NY,NX)=plt_rbgc%RUCNHB(N,L,NZ)
        RUNNOP(N,L,NZ,NY,NX)=plt_rbgc%RUNNOP(N,L,NZ)
        RUPNO3(N,L,NZ,NY,NX)=plt_rbgc%RUPNO3(N,L,NZ)
        RUONO3(N,L,NZ,NY,NX)=plt_rbgc%RUONO3(N,L,NZ)
        RUCNO3(N,L,NZ,NY,NX)=plt_rbgc%RUCNO3(N,L,NZ)
        RUNNXP(N,L,NZ,NY,NX)=plt_rbgc%RUNNXP(N,L,NZ)
        RUPNOB(N,L,NZ,NY,NX)=plt_rbgc%RUPNOB(N,L,NZ)
        RUONOB(N,L,NZ,NY,NX)=plt_rbgc%RUONOB(N,L,NZ)
        RUCNOB(N,L,NZ,NY,NX)=plt_rbgc%RUCNOB(N,L,NZ)
        RUPP2P(N,L,NZ,NY,NX)=plt_rbgc%RUPP2P(N,L,NZ)
        RUPH2P(N,L,NZ,NY,NX)=plt_rbgc%RUPH2P(N,L,NZ)
        RUOH2P(N,L,NZ,NY,NX)=plt_rbgc%RUOH2P(N,L,NZ)
        RUCH2P(N,L,NZ,NY,NX)=plt_rbgc%RUCH2P(N,L,NZ)
        RUPP2B(N,L,NZ,NY,NX)=plt_rbgc%RUPP2B(N,L,NZ)
        RUPH2B(N,L,NZ,NY,NX)=plt_rbgc%RUPH2B(N,L,NZ)
        RUOH2B(N,L,NZ,NY,NX)=plt_rbgc%RUOH2B(N,L,NZ)
        RUCH2B(N,L,NZ,NY,NX)=plt_rbgc%RUCH2B(N,L,NZ)
        RUPP1P(N,L,NZ,NY,NX)=plt_rbgc%RUPP1P(N,L,NZ)
        RUPH1P(N,L,NZ,NY,NX)=plt_rbgc%RUPH1P(N,L,NZ)
        RUOH1P(N,L,NZ,NY,NX)=plt_rbgc%RUOH1P(N,L,NZ)
        RUCH1P(N,L,NZ,NY,NX)=plt_rbgc%RUCH1P(N,L,NZ)
        RUPP1B(N,L,NZ,NY,NX)=plt_rbgc%RUPP1B(N,L,NZ)
        RUPH1B(N,L,NZ,NY,NX)=plt_rbgc%RUPH1B(N,L,NZ)
        RUOH1B(N,L,NZ,NY,NX)=plt_rbgc%RUOH1B(N,L,NZ)
        RUCH1B(N,L,NZ,NY,NX)=plt_rbgc%RUCH1B(N,L,NZ)
        ROXYP(N,L,NZ,NY,NX) =plt_rbgc%ROXYP(N,L,NZ)
        AllPlantRootH2OUptake_vr(N,L,NZ,NY,NX) =plt_ew%AllPlantRootH2OUptake_vr(N,L,NZ)
        WTRTL(N,L,NZ,NY,NX) =plt_biom%WTRTL(N,L,NZ)
        PopPlantRootC_vr(N,L,NZ,NY,NX) =plt_biom%PopPlantRootC_vr(N,L,NZ)
        WSRTL(N,L,NZ,NY,NX) =plt_biom%WSRTL(N,L,NZ)
        WFR(N,L,NZ,NY,NX)   =plt_rbgc%WFR(N,L,NZ)
      ENDDO
    ENDDO

    DO NR=1,JC
      NINR(NR,NZ,NY,NX)=plt_morph%NINR(NR,NZ)
      DO N=1,jroots
        RTWT1E(1:NumOfPlantChemElements,N,NR,NZ,NY,NX) =plt_biom%RTWT1E(1:NumOfPlantChemElements,N,NR,NZ)
        PrimRootDepth(N,NR,NZ,NY,NX) =plt_morph%PrimRootDepth(N,NR,NZ)
      ENDDO
      DO L=1,JZ
        DO N=1,2
          WTRT1E(1:NumOfPlantChemElements,N,L,NR,NZ,NY,NX) =plt_biom%WTRT1E(1:NumOfPlantChemElements,N,L,NR,NZ)
          WTRT2E(1:NumOfPlantChemElements,N,L,NR,NZ,NY,NX) =plt_biom%WTRT2E(1:NumOfPlantChemElements,N,L,NR,NZ)
          PrimRootLen(N,L,NR,NZ,NY,NX) =plt_morph%PrimRootLen(N,L,NR,NZ)
          SecndRootLen(N,L,NR,NZ,NY,NX) =plt_morph%SecndRootLen(N,L,NR,NZ)
          RTN2(N,L,NR,NZ,NY,NX)  =plt_morph%RTN2(N,L,NR,NZ)
        ENDDO
      ENDDO
    ENDDO

    EHVST(1:2,1:4,NZ,I,NY,NX)=plt_distb%EHVST(1:2,1:4,NZ)


    DO L=NU(NY,NX),NI(NZ,NY,NX)

      DO K=1,jcplx
        DO N=1,2
          DO NE=1,NumOfPlantChemElements
            RDFOME(NE,N,K,L,NZ,NY,NX)=plt_rbgc%RDFOME(NE,N,K,L,NZ)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO NE=1,NumOfPlantChemElements
      DO M=1,jsken
        DO N=0,JP
          CFOPE(NE,N,M,NZ,NY,NX)=plt_soilchem%CFOPE(NE,N,M,NZ)
        enddo
      enddo
    ENDDO
    MaxPrimRootRadius(2,NZ,NY,NX)=plt_morph%MaxPrimRootRadius(2,NZ)
    MaxSecndRootRadius(2,NZ,NY,NX)=plt_morph%MaxSecndRootRadius(2,NZ)
    RootPorosity(2,NZ,NY,NX)  =plt_morph%RootPorosity(2,NZ)
    UPMXZH(2,NZ,NY,NX)=plt_rbgc%UPMXZH(2,NZ)
    UPKMZH(2,NZ,NY,NX)=plt_rbgc%UPKMZH(2,NZ)
    UPMNZH(2,NZ,NY,NX)=plt_rbgc%UPMNZH(2,NZ)
    UPMXZO(2,NZ,NY,NX)=plt_rbgc%UPMXZO(2,NZ)
    UPKMZO(2,NZ,NY,NX)=plt_rbgc%UPKMZO(2,NZ)
    UPMNZO(2,NZ,NY,NX)=plt_rbgc%UPMNZO(2,NZ)
    UPMXPO(2,NZ,NY,NX)=plt_rbgc%UPMXPO(2,NZ)
    UPKMPO(2,NZ,NY,NX)=plt_rbgc%UPKMPO(2,NZ)
    UPMNPO(2,NZ,NY,NX)=plt_rbgc%UPMNPO(2,NZ)
    RSRR(2,NZ,NY,NX)  =plt_morph%RSRR(2,NZ)
    RSRA(2,NZ,NY,NX)  =plt_morph%RSRA(2,NZ)
    DO N=1,2
      RootPoreTortu4Gas(N,NZ,NY,NX)=plt_morph%RootPoreTortu4Gas(N,NZ)
      RRADP(N,NZ,NY,NX)=plt_morph%RRADP(N,NZ)
      DMVL(N,NZ,NY,NX) =plt_morph%DMVL(N,NZ)
      PrimRootSpecLen(N,NZ,NY,NX)=plt_morph%PrimRootSpecLen(N,NZ)
      SecndRootSpecLen(N,NZ,NY,NX)=plt_morph%SecndRootSpecLen(N,NZ)
      MaxPrimRootRadius1(N,NZ,NY,NX)=plt_morph%MaxPrimRootRadius1(N,NZ)
      MaxSecndRootRadius1(N,NZ,NY,NX)=plt_morph%MaxSecndRootRadius1(N,NZ)
      PrimRootXSecArea(N,NZ,NY,NX)=plt_morph%PrimRootXSecArea(N,NZ)
      SecndRootXSecArea(N,NZ,NY,NX)=plt_morph%SecndRootXSecArea(N,NZ)
    enDDO
  ENDDO

  end subroutine PlantAPIRecv


!------------------------------------------------------------------------------------------

  subroutine PlantAPISend(I,J,NY,NX)
  !
  !DESCRIPTION
  !Send data to plant model
  use PlantAPIData, only : plt_rad
  use EcoSIMConfig, only : jsken=>jskenc,jcplx => jcplxc
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer :: K,L,M,N,NB,NZ,NR,I1,NE

  plt_site%LYRC=LYRC
  I1=I+1;if(I1>LYRC)I1=1
  plt_site%ZERO=ZERO
  plt_site%ZERO2=ZERO2
  plt_site%ALAT=ALAT(NY,NX)
  plt_site%ATCA=ATCA(NY,NX)
  plt_morph%CanopyArea_grid=CanopyArea_grid(NY,NX)
  plt_morph%CanopyLA_grd=CanopyLA_grd(NY,NX)
  plt_site%ALT=ALT(NY,NX)
  plt_site%CCO2EI=CCO2EI(NY,NX)
  plt_site%CO2EI=CO2EI(NY,NX)
  plt_bgcr%CNETX=CNETX(NY,NX)
  plt_site%CO2E=CO2E(NY,NX)
  plt_site%COXYE=AtmGgms(idg_O2,NY,NX)
  plt_site%CCO2E=AtmGgms(idg_CO2,NY,NX)
  plt_site%CCH4E=AtmGgms(idg_CH4,NY,NX)
  plt_site%CZ2OE=AtmGgms(idg_N2O,NY,NX)
  plt_site%CH2GE=AtmGgms(idg_H2,NY,NX)
  plt_site%CNH3E=AtmGgms(idg_NH3,NY,NX)
  plt_site%DYLX=DYLX(NY,NX)
  plt_site%DYLN=DYLN(NY,NX)
  plt_ew%SnowDepth=SnowDepth(NY,NX)
  plt_site%DYLM=DYLM(NY,NX)
  plt_site%IETYP=IETYP(NY,NX)
  plt_site%IYRC=IYRC
  plt_site%NL=NL(NY,NX)
  plt_site%NP0=NP0(NY,NX)
  plt_site%NJ=NJ(NY,NX)
  plt_site%NP=NP(NY,NX)
  plt_site%NU=NU(NY,NX)
  plt_site%OXYE=OXYE(NY,NX)
  plt_ew%BndlResistAboveCanG=BndlResistAboveCanG(NY,NX)
  plt_ew%RIB=RIB(NY,NX)
  plt_rad%SSINN=SSINN(NY,NX)
  plt_rad%SSIN=SSIN(NY,NX)
  plt_ew%TKSnow=TKSnow(1,NY,NX)  !surface layer snow temperature
  plt_ew%TairK=TairK(NY,NX)
  plt_rad%LWRadGrnd=LWRadGrnd(NY,NX)
  plt_rad%LWRadSky=LWRadSky(NY,NX)
  plt_ew%VPA=VPA(NY,NX)
  plt_distb%XCORP=XCORP(NY,NX)
  plt_site%ZNOON =ZNOON(NY,NX)
  plt_site%ZEROS2=ZEROS2(NY,NX)
  plt_site%ZEROS =ZEROS(NY,NX)
  plt_ew%RoughHeight=RoughHeight(NY,NX)
  plt_morph%GridMaxCanopyHeight=GridMaxCanopyHeight(NY,NX)
  plt_ew%ZeroPlanDisp=ZeroPlanDisp(NY,NX)
  plt_site%IDATA(:)=IDATA(:)
  plt_distb%DCORP=DCORP(I,NY,NX)
  plt_distb%ITILL=ITILL(I,NY,NX)

  plt_morph%CanopyHeightz(0)=CanopyHeightz(0,NY,NX)
  DO  L=1,JC
    plt_morph%CanopyHeightz(L)=CanopyHeightz(L,NY,NX)
    plt_rad%TAUS(L)=TAUS(L,NY,NX)
    plt_rad%TAU0(L)=TAU0(L,NY,NX)
  ENDDO
  plt_rad%TAUS(JC+1)=TAUS(JC+1,NY,NX)
  plt_rad%TAU0(JC+1)=TAU0(JC+1,NY,NX)

  DO N=1,JLI
    plt_rad%ZSIN(N)=ZSIN(N)
  ENDDO

  DO L=1,NL(NY,NX)
    plt_soilchem%HydroCondMicP4RootUptake(L)=HydroCondMicP4RootUptake(L,NY,NX)
    plt_soilchem%GasDifc(idg_beg:idg_end,L)=GasDifc(idg_beg:idg_end,L,NY,NX)
    plt_soilchem%SoilResit4RootPentration(L)=SoilResit4RootPentration(L,NY,NX)
    plt_site%DPTHZ(L)=DPTHZ(L,NY,NX)
  ENDDO
  plt_allom%FVRN(:) = FVRN(:)
  DO L=1,NL(NY,NX)
    plt_soilchem%trc_gasml(idg_CO2,L)=trc_gasml(idg_CO2,L,NY,NX)
    plt_soilchem%trc_gasml(idg_O2,L)=trc_gasml(idg_O2,L,NY,NX)
  ENDDO

  DO L=0,NL(NY,NX)
    plt_site%AREA3(L)     =AREA(3,L,NY,NX)
    plt_soilchem%SoiBulkDensity(L)  =SoiBulkDensity(L,NY,NX)

    plt_soilchem%trc_solcl(ids_beg:ids_end,L) =trc_solcl(ids_beg:ids_end,L,NY,NX)

    plt_soilchem%SolDifc(ids_beg:ids_end,L) =SolDifc(ids_beg:ids_end,L,NY,NX)

    plt_soilchem%trc_gascl(idg_beg:idg_end,L) =trc_gascl(idg_beg:idg_end,L,NY,NX)
    plt_soilchem%CORGC(L) =CORGC(L,NY,NX)
    plt_site%CumSoilThickness(L)    =CumSoilThickness(L,NY,NX)
    plt_site%FracSoiAsMicP(L) =FracSoiAsMicP(L,NY,NX)

    plt_soilchem%trc_solml(ids_beg:ids_end,L)  =trc_solml(ids_beg:ids_end,L,NY,NX)

    plt_soilchem%GasSolbility(idg_beg:idg_end-1,L) =GasSolbility(idg_beg:idg_end-1,L,NY,NX)

    plt_ew%TotalSoilH2OPSIMPa(L)       =TotalSoilH2OPSIMPa(L,NY,NX)
    plt_bgcr%RPO4Y(L)     =RPO4Y(L,NY,NX)
    plt_bgcr%RPOBY(L)     =RPOBY(L,NY,NX)
    plt_bgcr%RP14Y(L)     =RP14Y(L,NY,NX)
    plt_bgcr%RP1BY(L)     =RP1BY(L,NY,NX)
    plt_bgcr%RNO3Y(L)     =RNO3Y(L,NY,NX)
    plt_bgcr%RNH4Y(L)     =RNH4Y(L,NY,NX)
    plt_bgcr%RNHBY(L)     =RNHBY(L,NY,NX)
    plt_bgcr%RN3BY(L)     =RN3BY(L,NY,NX)
    plt_bgcr%ROXYF(L)     =ROXYF(L,NY,NX)
    plt_bgcr%RCO2F(L)     =RCO2F(L,NY,NX)
    plt_bgcr%ROXYL(L)     =ROXYL(L,NY,NX)
    plt_bgcr%ROXYY(L)     =ROXYY(L,NY,NX)
    plt_ew%TKS(L)         =TKS(L,NY,NX)


    plt_soilchem%THETW(L) =THETW(L,NY,NX)
    plt_soilchem%THETY(L) =THETY(L,NY,NX)
    plt_soilchem%TFND(L)  =TFND(L,NY,NX)
    plt_soilchem%VLSoilPoreMicP(L)  =VLSoilPoreMicP(L,NY,NX)
    plt_soilchem%trcs_VLN(ids_H1PO4B,L) =trcs_VLN(ids_H1PO4B,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NO3,L) =trcs_VLN(ids_NO3,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_H1PO4,L) =trcs_VLN(ids_H1PO4,L,NY,NX)
    plt_soilchem%VLSoilMicP(L)  =VLSoilMicP(L,NY,NX)
    plt_soilchem%VLiceMicP(L)  =VLiceMicP(L,NY,NX)
    plt_soilchem%VLWatMicP(L)  =VLWatMicP(L,NY,NX)
    plt_soilchem%VLMicP(L)  =VLMicP(L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NO3B,L) =trcs_VLN(ids_NO3B,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NH4,L) =trcs_VLN(ids_NH4,L,NY,NX)
    plt_soilchem%trcs_VLN(ids_NH4B,L) =trcs_VLN(ids_NH4B,L,NY,NX)

    plt_site%DLYR3(L)     =DLYR(3,L,NY,NX)
    DO K=1,jcplx
      plt_soilchem%FOSRH(K,L)=FOSRH(K,L,NY,NX)
      plt_soilchem%DOM(idom_doc,K,L)=DOM(idom_doc,K,L,NY,NX)
      plt_soilchem%DOM(idom_don,K,L)=DOM(idom_don,K,L,NY,NX)
      plt_soilchem%DOM(idom_dop,K,L)=DOM(idom_dop,K,L,NY,NX)
    ENDDO
  ENDDO

  DO NZ=1,NP0(NY,NX)
!plant properties begin
    plt_photo%ICTYP(NZ)=ICTYP(NZ,NY,NX)
    plt_pheno%IGTYP(NZ)=IGTYP(NZ,NY,NX)
    plt_pheno%ISTYP(NZ)=ISTYP(NZ,NY,NX)
    plt_pheno%IDTYP(NZ)=IDTYP(NZ,NY,NX)
    plt_pheno%IWTYP(NZ)=IWTYP(NZ,NY,NX)
    plt_pheno%IPTYP(NZ)=IPTYP(NZ,NY,NX)
    plt_pheno%IBTYP(NZ)=IBTYP(NZ,NY,NX)
    plt_pheno%ZTYPI(NZ)=ZTYPI(NZ,NY,NX)
    plt_morph%IRTYP(NZ)=IRTYP(NZ,NY,NX)
    plt_morph%INTYP(NZ)=INTYP(NZ,NY,NX)
    plt_morph%MY(NZ)=MY(NZ,NY,NX)
    plt_photo%VCMX(NZ)=VCMX(NZ,NY,NX)
    plt_photo%VOMX(NZ)=VOMX(NZ,NY,NX)
    plt_photo%VCMX4(NZ)=VCMX4(NZ,NY,NX)
    plt_photo%XKCO2(NZ)=XKCO2(NZ,NY,NX)
    plt_photo%XKO2(NZ)=XKO2(NZ,NY,NX)
    plt_photo%XKCO24(NZ)=XKCO24(NZ,NY,NX)
    plt_photo%RUBP(NZ)=RUBP(NZ,NY,NX)
    plt_photo%PEPC(NZ)=PEPC(NZ,NY,NX)
    plt_photo%ETMX(NZ)=ETMX(NZ,NY,NX)
    plt_photo%CHL(NZ)=CHL(NZ,NY,NX)
    plt_photo%CHL4(NZ)=CHL4(NZ,NY,NX)
    plt_photo%CanPCi2CaRatio(NZ)=CanPCi2CaRatio(NZ,NY,NX)

    plt_pheno%XRNI(NZ)=XRNI(NZ,NY,NX)
    plt_pheno%XRLA(NZ)=XRLA(NZ,NY,NX)
    plt_pheno%CTC(NZ)=CTC(NZ,NY,NX)
    plt_morph%WDLF(NZ)=WDLF(NZ,NY,NX)
    plt_pheno%PB(NZ)=PB(NZ,NY,NX)
    plt_morph%XTLI(NZ)=XTLI(NZ,NY,NX)
    plt_pheno%XDL(NZ)=XDL(NZ,NY,NX)
    plt_pheno%XPPD(NZ)=XPPD(NZ,NY,NX)
    plt_morph%SLA1(NZ)=SLA1(NZ,NY,NX)
    plt_morph%SSL1(NZ)=SSL1(NZ,NY,NX)
    plt_morph%SNL1(NZ)=SNL1(NZ,NY,NX)
    DO  N=1,JLI
      plt_morph%CLASS(N,NZ)=CLASS(N,NZ,NY,NX)
    ENDDO
    plt_morph%ClumpFactort0(NZ)=ClumpFactort0(NZ,NY,NX)
    plt_morph%ANGBR(NZ)=ANGBR(NZ,NY,NX)
    plt_morph%ANGSH(NZ)=ANGSH(NZ,NY,NX)
    plt_morph%STMX(NZ)=STMX(NZ,NY,NX)
    plt_morph%SDMX(NZ)=SDMX(NZ,NY,NX)
    plt_morph%MaxSeedCMass(NZ)=MaxSeedCMass(NZ,NY,NX)
    plt_morph%SeedCMass(NZ)=SeedCMass(NZ,NY,NX)
    plt_pheno%GFILL(NZ)=GFILL(NZ,NY,NX)
    plt_biom%WTSTDI(NZ)=WTSTDI(NZ,NY,NX)

!initial root values
    DO N=1,MY(NZ,NY,NX)
      plt_morph%MaxPrimRootRadius(N,NZ)=MaxPrimRootRadius(N,NZ,NY,NX)
      plt_morph%MaxSecndRootRadius(N,NZ)=MaxSecndRootRadius(N,NZ,NY,NX)
      plt_morph%RootPorosity(N,NZ)  =RootPorosity(N,NZ,NY,NX)
      plt_morph%RSRR(N,NZ)  =RSRR(N,NZ,NY,NX)
      plt_morph%RSRA(N,NZ)=RSRA(N,NZ,NY,NX)
      plt_rbgc%UPMXZH(N,NZ)=UPMXZH(N,NZ,NY,NX)
      plt_rbgc%UPKMZH(N,NZ)=UPKMZH(N,NZ,NY,NX)
      plt_rbgc%UPMNZH(N,NZ)=UPMNZH(N,NZ,NY,NX)
      plt_rbgc%UPMXZO(N,NZ)=UPMXZO(N,NZ,NY,NX)
      plt_rbgc%UPKMZO(N,NZ)=UPKMZO(N,NZ,NY,NX)
      plt_rbgc%UPMNZO(N,NZ)=UPMNZO(N,NZ,NY,NX)
      plt_rbgc%UPMXPO(N,NZ)=UPMXPO(N,NZ,NY,NX)
      plt_rbgc%UPKMPO(N,NZ)=UPKMPO(N,NZ,NY,NX)
      plt_rbgc%UPMNPO(N,NZ)=UPMNPO(N,NZ,NY,NX)
    ENDDO
    plt_pheno%PR(NZ)=PR(NZ,NY,NX)
    plt_pheno%PTSHT(NZ)=PTSHT(NZ,NY,NX)
    plt_morph%RTFQ(NZ)=RTFQ(NZ,NY,NX)
    plt_ew%OSMO(NZ)=OSMO(NZ,NY,NX)
    plt_photo%RCS(NZ)=RCS(NZ,NY,NX)
    plt_photo%RSMX(NZ)=RSMX(NZ,NY,NX)

    plt_allom%DMLF(NZ) =DMLF(NZ,NY,NX)
    plt_allom%DMSHE(NZ)=DMSHE(NZ,NY,NX)
    plt_allom%DMSTK(NZ)=DMSTK(NZ,NY,NX)
    plt_allom%DMRSV(NZ)=DMRSV(NZ,NY,NX)
    plt_allom%DMHSK(NZ)=DMHSK(NZ,NY,NX)
    plt_allom%DMEAR(NZ)=DMEAR(NZ,NY,NX)
    plt_allom%DMGR(NZ) =DMGR(NZ,NY,NX)
    plt_allom%BiomGrowthYieldRoot(NZ) =BiomGrowthYieldRoot(NZ,NY,NX)
    plt_allom%DMND(NZ) =DMND(NZ,NY,NX)
    plt_allom%CNLF(NZ) =CNLF(NZ,NY,NX)
    plt_allom%CNSHE(NZ)=CNSHE(NZ,NY,NX)
    plt_allom%CNSTK(NZ)=CNSTK(NZ,NY,NX)
    plt_allom%CNRSV(NZ)=CNRSV(NZ,NY,NX)
    plt_allom%CNHSK(NZ)=CNHSK(NZ,NY,NX)
    plt_allom%CNEAR(NZ)=CNEAR(NZ,NY,NX)
    plt_allom%CNGR(NZ)=CNGR(NZ,NY,NX)
    plt_allom%CNRT(NZ)=CNRT(NZ,NY,NX)
    plt_allom%CNND(NZ)=CNND(NZ,NY,NX)
    plt_allom%CPLF(NZ)=CPLF(NZ,NY,NX)
    plt_allom%CPSHE(NZ)=CPSHE(NZ,NY,NX)
    plt_allom%CPSTK(NZ)=CPSTK(NZ,NY,NX)
    plt_allom%CPRSV(NZ)=CPRSV(NZ,NY,NX)
    plt_allom%CPHSK(NZ)=CPHSK(NZ,NY,NX)
    plt_allom%CPEAR(NZ)=CPEAR(NZ,NY,NX)
    plt_allom%CPGR(NZ)=CPGR(NZ,NY,NX)
    plt_allom%CPRT(NZ)=CPRT(NZ,NY,NX)
    plt_allom%CPND(NZ)=CPND(NZ,NY,NX)

!plant properties end

    plt_morph%CanopyArea_pft(NZ)=CanopyArea_pft(NZ,NY,NX)
    plt_distb%IYRX(NZ)=IYRX(NZ,NY,NX)
    plt_distb%IDAYX(NZ)=IDAYX(NZ,NY,NX)
    plt_distb%IYRY(NZ)=IYRY(NZ,NY,NX)
    plt_rad%PARByCanP(NZ)=PARByCanP(NZ,NY,NX)
    plt_rad%SWRadByCanP(NZ)=SWRadByCanP(NZ,NY,NX)
    plt_ew%PrecIntcptByCanP(NZ)=PrecIntcptByCanP(NZ,NY,NX)

    plt_site%PPZ(NZ)=PPZ(NZ,NY,NX)
    plt_distb%IDAYY(NZ)=IDAYY(NZ,NY,NX)
    plt_morph%ClumpFactort(NZ)=ClumpFactort(NZ,NY,NX)
    plt_site%DATAP(NZ)=DATAP(NZ,NY,NX)
    plt_pheno%GROUPI(NZ)=GROUPI(NZ,NY,NX)
    plt_biom%WTSHTA(NZ)=WTSHTA(NZ,NY,NX)


    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      plt_pheno%VRNL(NB,NZ)=VRNL(NB,NZ,NY,NX)
      plt_pheno%VRNX(NB,NZ)=VRNX(NB,NZ,NY,NX)
    ENDDO

    DO L=1,JC
      DO  M=1,JSA
        DO  N=1,JLI
          plt_rad%PAR(N,M,L,NZ)=PAR(N,M,L,NZ,NY,NX)
          plt_rad%PARDIF(N,M,L,NZ)=PARDIF(N,M,L,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,NL(NY,NX)
    DO M=1,NPH
      plt_site%VLWatMicPM(M,L)=VLWatMicPM(M,L,NY,NX)
      plt_site%VLsoiAirPM(M,L)=VLsoiAirPM(M,L,NY,NX)
      plt_site%TortMicPM(M,L)=TortMicPM(M,L,NY,NX)
      plt_site%FILM(M,L)=FILM(M,L,NY,NX)
      plt_soilchem%DiffusivitySolutEff(M,L)=DiffusivitySolutEff(M,L,NY,NX)
    ENDDO
  ENDDO

! sent variables also modified
  plt_site%IFLGT=IFLGT(NY,NX)
  plt_site%VOLWOU=VOLWOU
  plt_site%PPT=PPT(NY,NX)
  plt_bgcr%RECO=RECO(NY,NX)
  plt_biom%WTSTGET(1:NumOfPlantChemElements)=WTSTGET(1:NumOfPlantChemElements,NY,NX)
  plt_bgcr%ZESNC(1:NumOfPlantChemElements)=ZESNC(1:NumOfPlantChemElements,NY,NX)
  plt_morph%StemAreag=StemAreag(NY,NX)
  plt_ew%TSH=TSH(NY,NX)
  plt_ew%CanH2OHeldVg=CanH2OHeldVg(NY,NX)
  plt_bgcr%Eco_NBP_col=Eco_NBP_col(NY,NX)
  plt_bgcr%Eco_GPP_col=Eco_GPP_col(NY,NX)
  plt_ew%TLEC=TLEC(NY,NX)
  plt_ew%TSHC=TSHC(NY,NX)
  plt_bgcr%Eco_AutoR_col=Eco_AutoR_col(NY,NX)
  plt_ew%UVOLO    =UVOLO(NY,NX)
  plt_distb%XHVSTE(1:NumOfPlantChemElements)=XHVSTE(1:NumOfPlantChemElements,NY,NX)
  plt_rad%TRN     =TRN(NY,NX)
  plt_ew%VapXAir2CanG=VapXAir2CanG(NY,NX)
  plt_ew%TLE=TLE(NY,NX)
  plt_rbgc%TRootGasLoss_disturb(idg_beg:idg_end-1)=TRootGasLoss_disturb(idg_beg:idg_end-1,NY,NX)
  plt_ew%TGH    =TGH(NY,NX)
  plt_ew%TEVAPP =TEVAPP(NY,NX)
  plt_ew%THFLXC =THFLXC(NY,NX)
  plt_ew%LWRadCanG  =LWRadCanG(NY,NX)
  plt_ew%CanWatg =CanWatg(NY,NX)
  plt_ew%TENGYC =TENGYC(NY,NX)
  plt_bgcr%TCCAN=TCCAN(NY,NX)
  plt_distb%FERT(1:20)=FERT(1:20,I1,NY,NX)

  DO  L=1,JC
    plt_morph%CanopyStemA_lyr(L)=CanopyStemA_lyr(L,NY,NX)
    plt_biom%WGLFT(L)=WGLFT(L,NY,NX)
    plt_morph%CanopyLAgrid_lyr(L)=CanopyLAgrid_lyr(L,NY,NX)
  ENDDO

  plt_allom%FWOODE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)=FWOODE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)
  plt_allom%FWODRE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)=FWODRE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)
  plt_allom%FWODBE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)=FWODBE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)
  plt_allom%FWODLE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)=FWODLE(1:NumOfPlantChemElements,1:NumOfPlantLitrCmplxs)

  DO L=0,NL(NY,NX)
    DO K=1,jcplx
      plt_bgcr%RDOM_micb_flx(idom_doc,K,L)=RDOM_micb_flx(idom_doc,K,L,NY,NX)
      plt_bgcr%RDOM_micb_flx(idom_don,K,L)=RDOM_micb_flx(idom_don,K,L,NY,NX)
      plt_bgcr%RDOM_micb_flx(idom_dop,K,L)=RDOM_micb_flx(idom_dop,K,L,NY,NX)
    ENDDO
  ENDDO

  DO L=0,JZ
    plt_bgcr%RPOBX(L)=RPOBX(L,NY,NX)
    plt_bgcr%RP1BX(L)=RP1BX(L,NY,NX)
    plt_bgcr%RN3BX(L)=RN3BX(L,NY,NX)
    plt_bgcr%RNHBX(L)=RNHBX(L,NY,NX)
    plt_bgcr%RP14X(L)=RP14X(L,NY,NX)
    plt_bgcr%RPO4X(L)=RPO4X(L,NY,NX)
    plt_bgcr%RNO3X(L)=RNO3X(L,NY,NX)
    plt_bgcr%RNH4X(L)=RNH4X(L,NY,NX)
    plt_bgcr%ROXYX(L)=ROXYX(L,NY,NX)
    plt_ew%THeatRootUptake(L)  =THeatRootUptake(L,NY,NX)
    plt_ew%GridPlantRootH2OUptake_vr(L) =GridPlantRootH2OUptake_vr(L,NY,NX)
    DO  K=1,micpar%NumOfPlantLitrCmplxs
      DO NE=1,NumOfPlantChemElements
        DO  M=1,jsken
          plt_bgcr%LitrfalChemElemnts_vr(NE,M,K,L)=LitrfalChemElemnts_vr(NE,M,K,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO L=1,JZ
    plt_rbgc%TUPNF(L)=TUPNF(L,NY,NX)
    plt_morph%RTDNT(L)=RTDNT(L,NY,NX)
    plt_rbgc%trcg_TLP(idg_beg:idg_end-1,L)=trcg_TLP(idg_beg:idg_end-1,L,NY,NX)
    plt_rbgc%trcg_TFLA(idg_beg:idg_end-1,L)=trcg_TFLA(idg_beg:idg_end-1,L,NY,NX)
    plt_bgcr%TCO2P(L) =TCO2P(L,NY,NX)
    plt_bgcr%TUPOXP(L)=TUPOXP(L,NY,NX)
    plt_bgcr%TCO2S(L) =TCO2S(L,NY,NX)
    plt_bgcr%TUPOXS(L)=TUPOXS(L,NY,NX)
    plt_bgcr%TUPCHS(L)=TUPCHS(L,NY,NX)
    plt_bgcr%TUPN2S(L)=TUPN2S(L,NY,NX)
    plt_bgcr%TUPN3S(L)=TUPN3S(L,NY,NX)
    plt_bgcr%TUPN3B(L)=TUPN3B(L,NY,NX)
    plt_bgcr%TUPHGS(L)=TUPHGS(L,NY,NX)
    plt_bgcr%TUPNH4(L)=TUPNH4(L,NY,NX)
    plt_bgcr%TUPNO3(L)=TUPNO3(L,NY,NX)
    plt_bgcr%TUPH1B(L)=TUPH1B(L,NY,NX)
    plt_bgcr%TUPH2B(L)=TUPH2B(L,NY,NX)
    plt_bgcr%TUPNHB(L)=TUPNHB(L,NY,NX)
    plt_bgcr%TUPNOB(L)=TUPNOB(L,NY,NX)
    plt_bgcr%TUPH1P(L)=TUPH1P(L,NY,NX)
    plt_bgcr%TUPH2P(L)=TUPH2P(L,NY,NX)
    DO  K=1,jcplx
      plt_bgcr%TDFOME(1:NumOfPlantChemElements,K,L)=TDFOME(1:NumOfPlantChemElements,K,L,NY,NX)
    ENDDO
  ENDDO


  DO NZ=1,NP0(NY,NX)
    plt_biom%WTRTE(1:NumOfPlantChemElements,NZ)=WTRTE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTLFE(1:NumOfPlantChemElements,NZ)  =WTLFE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTEARE(1:NumOfPlantChemElements,NZ) =WTEARE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTRTSE(1:NumOfPlantChemElements,NZ) =WTRTSE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTHSKE(1:NumOfPlantChemElements,NZ) =WTHSKE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTRSVE(1:NumOfPlantChemElements,NZ) =WTRSVE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTSTKE(1:NumOfPlantChemElements,NZ) =WTSTKE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTSHEE(1:NumOfPlantChemElements,NZ) =WTSHEE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTGRE(1:NumOfPlantChemElements,NZ)=WTGRE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_site%BALE(1:NumOfPlantChemElements,NZ)  =BALE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ)=CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%CanopyNonstructElements_pft(1:NumOfPlantChemElements,NZ)=CanopyNonstructElements_pft(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%EPOLNP(1:NumOfPlantChemElements,NZ)=EPOLNP(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_bgcr%HESNC(1:NumOfPlantChemElements,NZ) =HESNC(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_distb%HVSTE(1:NumOfPlantChemElements,NZ)=HVSTE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_pheno%RSETE(1:NumOfPlantChemElements,NZ)=RSETE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_bgcr%TESN0(1:NumOfPlantChemElements,NZ)  =TESN0(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_bgcr%TESNC(1:NumOfPlantChemElements,NZ)  =TESNC(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_rbgc%TEUPTK(1:NumOfPlantChemElements,NZ) =TEUPTK(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_distb%THVSTE(1:NumOfPlantChemElements,NZ)=THVSTE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_rbgc%UPOME(1:NumOfPlantChemElements,NZ)=UPOME(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTRVE(1:NumOfPlantChemElements,NZ) =WTRVE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%CanPShootElmMass(1:NumOfPlantChemElements,NZ) =CanPShootElmMass(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTSTGE(1:NumOfPlantChemElements,NZ) =WTSTGE(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_biom%WTNDE(1:NumOfPlantChemElements,NZ)  =WTNDE(1:NumOfPlantChemElements,NZ,NY,NX)

    plt_ew%TKCZ(NZ)=TKCZ(NZ,NY,NX)
    plt_photo%SO2(NZ)=SO2(NZ,NY,NX)
    plt_ew%PSICanPDailyMin(NZ)=PSICanPDailyMin(NZ,NY,NX)
    plt_ew%HeatStorCanP(NZ)=HeatStorCanP(NZ,NY,NX)
    plt_photo%CO2Q(NZ)=CO2Q(NZ,NY,NX)
    plt_photo%O2L(NZ)=O2L(NZ,NY,NX)

    plt_photo%CO2L(NZ)=CO2L(NZ,NY,NX)
    plt_distb%EHVST(1:2,1:4,NZ)=EHVST(1:2,1:4,NZ,I,NY,NX)

    plt_pheno%IDTH(NZ)=IDTH(NZ,NY,NX)
    plt_distb%IYR0(NZ)=IYR0(NZ,NY,NX)
    plt_morph%NNOD(NZ)=NNOD(NZ,NY,NX)

    plt_pheno%IDTHP(NZ)=IDTHP(NZ,NY,NX)
    plt_pheno%IDTHR(NZ)=IDTHR(NZ,NY,NX)
    plt_allom%CNRTS(NZ)=CNRTS(NZ,NY,NX)
    plt_allom%CPRTS(NZ)=CPRTS(NZ,NY,NX)

    plt_rbgc%ZEROQ(NZ)=ZEROQ(NZ,NY,NX)
    plt_pheno%SSTX(NZ)=SSTX(NZ,NY,NX)
    plt_rad%FracPARByCanP(NZ) =FracPARByCanP(NZ,NY,NX)
    plt_bgcr%RNH3C(NZ)=RNH3C(NZ,NY,NX)
    plt_ew%DTKC(NZ)   =DTKC(NZ,NY,NX)
    plt_pheno%ZTYP(NZ)=ZTYP(NZ,NY,NX)
    plt_pheno%HTC(NZ) =HTC(NZ,NY,NX)
    plt_ew%PTrans(NZ)      =PTrans(NZ,NY,NX)

    plt_biom%CanPStalkC(NZ) =CanPStalkC(NZ,NY,NX)

    plt_photo%CHILL(NZ)=CHILL(NZ,NY,NX)
    plt_ew%PSICanPOsmo(NZ)=PSICanPOsmo(NZ,NY,NX)

    plt_ew%TCC(NZ)=TCC(NZ,NY,NX)
    plt_allom%RootFracRemobilizableBiom(NZ)=RootFracRemobilizableBiom(NZ,NY,NX)
    plt_photo%MaxCanPStomaResistH2O(NZ)=MaxCanPStomaResistH2O(NZ,NY,NX)
    plt_biom%WTRTA(NZ)=WTRTA(NZ,NY,NX)
    plt_morph%ClumpFactor(NZ)=ClumpFactor(NZ,NY,NX)

    plt_site%PPI(NZ)=PPI(NZ,NY,NX)
    plt_site%PPX(NZ)=PPX(NZ,NY,NX)
    plt_distb%VOXYF(NZ)=VOXYF(NZ,NY,NX)
    plt_ew%ENGYX(NZ)=ENGYX(NZ,NY,NX)
    plt_bgcr%TCO2A(NZ)=TCO2A(NZ,NY,NX)
    plt_ew%ETCanP(NZ)=ETCanP(NZ,NY,NX)

    plt_morph%NBT(NZ)=NBT(NZ,NY,NX)
    plt_morph%NGTopRootLayer(NZ)=NGTopRootLayer(NZ,NY,NX)
    plt_pheno%IFLGI(NZ)=IFLGI(NZ,NY,NX)
    plt_morph%NIXBotRootLayer(NZ)=NIXBotRootLayer(NZ,NY,NX)
    plt_morph%NRT(NZ)=NRT(NZ,NY,NX)
    plt_morph%NB1(NZ)=NB1(NZ,NY,NX)
    plt_morph%NumOfBranches_pft(NZ)=NumOfBranches_pft(NZ,NY,NX)

    plt_pheno%IFLGC(NZ)=IFLGC(NZ,NY,NX)
    plt_distb%IDAY0(NZ)=IDAY0(NZ,NY,NX)
    plt_distb%IDAYH(NZ)=IDAYH(NZ,NY,NX)
    plt_distb%IYRH(NZ)=IYRH(NZ,NY,NX)


    plt_distb%HVST(NZ) =HVST(NZ,I,NY,NX)
    plt_distb%IHVST(NZ)=IHVST(NZ,I,NY,NX)
    plt_distb%JHVST(NZ)=JHVST(NZ,I,NY,NX)
    plt_distb%THIN_pft(NZ) =THIN_pft(NZ,I,NY,NX)
    plt_morph%CanPSA(NZ)=CanPSA(NZ,NY,NX)
    plt_morph%CanopyLeafA_pft(NZ)=CanopyLeafA_pft(NZ,NY,NX)

    plt_photo%O2I(NZ)  =O2I(NZ,NY,NX)
    plt_photo%CO2I(NZ) =CO2I(NZ,NY,NX)

    plt_biom%NoduleNonstructCconc_pft(NZ)=NoduleNonstructCconc_pft(NZ,NY,NX)
    plt_bgcr%CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ,NY,NX)
    plt_bgcr%CARBN(NZ)=CARBN(NZ,NY,NX)

    plt_allom%CNWS(NZ)=CNWS(NZ,NY,NX)
    plt_allom%CPWS(NZ)=CPWS(NZ,NY,NX)

    plt_photo%DCO2(NZ)=DCO2(NZ,NY,NX)
    plt_photo%FMOL(NZ)=FMOL(NZ,NY,NX)
    plt_allom%FNOD(NZ)=FNOD(NZ,NY,NX)

    plt_morph%HypoctoylHeight(NZ)=HypoctoylHeight(NZ,NY,NX)
    plt_rbgc%HEUPTK(1:NumOfPlantChemElements,NZ)=HEUPTK(1:NumOfPlantChemElements,NZ,NY,NX)
    plt_morph%CanPHeight4WatUptake(NZ)=CanPHeight4WatUptake(NZ,NY,NX)
    plt_morph%NI(NZ)   =NI(NZ,NY,NX)
    plt_photo%CanPbndlResist(NZ)   =CanPbndlResist(NZ,NY,NX)
    plt_photo%CanPStomaResistH2O(NZ)   =CanPStomaResistH2O(NZ,NY,NX)
    plt_ew%TKC(NZ)     =TKC(NZ,NY,NX)
    plt_ew%HeatXAir2PCan(NZ)   =HeatXAir2PCan(NZ,NY,NX)
    plt_rad%RadNet2CanP(NZ)   =RadNet2CanP(NZ,NY,NX)
    plt_rad%LWRadCanP(NZ)  =LWRadCanP(NZ,NY,NX)
    plt_ew%EvapTransHeatP(NZ)   =EvapTransHeatP(NZ,NY,NX)
    plt_ew%VapXAir2PCan(NZ)   =VapXAir2PCan(NZ,NY,NX)
    plt_photo%MinCanPStomaResistH2O(NZ) =MinCanPStomaResistH2O(NZ,NY,NX)
    plt_pheno%OFFST(NZ)=OFFST(NZ,NY,NX)
    plt_pheno%PlantO2Stress(NZ)=PlantO2Stress(NZ,NY,NX)
    plt_site%pftPlantPopulation(NZ)=pftPlantPopulation(NZ,NY,NX)
    plt_ew%PSICanP(NZ)=PSICanP(NZ,NY,NX)
    plt_ew%PSICanPTurg(NZ)=PSICanPTurg(NZ,NY,NX)
    plt_photo%RCMX(NZ)=RCMX(NZ,NY,NX)
    plt_bgcr%RootGasLoss_disturb(idg_beg:idg_end-1,NZ)=RootGasLoss_disturb(idg_beg:idg_end-1,NZ,NY,NX)
    plt_ew%RAZ(NZ)=RAZ(NZ,NY,NX)
    plt_photo%SCO2(NZ)=SCO2(NZ,NY,NX)
    plt_morph%SeedinDepth(NZ)=SeedinDepth(NZ,NY,NX)
    plt_morph%PlantinDepth(NZ)=PlantinDepth(NZ,NY,NX)
    plt_morph%SeedLength(NZ)=SeedLength(NZ,NY,NX)
    plt_morph%SeedVolume(NZ)=SeedVolume(NZ,NY,NX)
    plt_morph%SeedArea(NZ)=SeedArea(NZ,NY,NX)
    plt_pheno%TCZ(NZ)    =TCZ(NZ,NY,NX)
    plt_pheno%TCG(NZ)    =TCG(NZ,NY,NX)
    plt_pheno%TCX(NZ)    =TCX(NZ,NY,NX)
    plt_pheno%TKG(NZ)      =TKG(NZ,NY,NX)
    plt_pheno%fTgrowCanP(NZ)  =fTgrowCanP(NZ,NY,NX)

    plt_bgcr%TCO2T(NZ)  =TCO2T(NZ,NY,NX)
    plt_photo%XKCO2L(NZ)=XKCO2L(NZ,NY,NX)
    plt_photo%XKCO2O(NZ)=XKCO2O(NZ,NY,NX)
    plt_bgcr%TNH3C(NZ)  =TNH3C(NZ,NY,NX)
    plt_bgcr%TZUPFX(NZ) =TZUPFX(NZ,NY,NX)
    plt_rbgc%UPNF(NZ)=UPNF(NZ,NY,NX)
    plt_rbgc%UPNO3(NZ)=UPNO3(NZ,NY,NX)
    plt_rbgc%UPNH4(NZ)=UPNH4(NZ,NY,NX)
    plt_rbgc%UPH1P(NZ)=UPH1P(NZ,NY,NX)
    plt_rbgc%UPH2P(NZ)=UPH2P(NZ,NY,NX)
    plt_ew%WatByPCan(NZ)=WatByPCan(NZ,NY,NX)
    plt_ew%VHeatCapCanP(NZ)=VHeatCapCanP(NZ,NY,NX)
    plt_distb%VCH4F(NZ)=VCH4F(NZ,NY,NX)
    plt_distb%VCO2F(NZ)=VCO2F(NZ,NY,NX)
    plt_distb%VN2OF(NZ)=VN2OF(NZ,NY,NX)
    plt_distb%VNH3F(NZ)=VNH3F(NZ,NY,NX)
    plt_distb%VPO4F(NZ)=VPO4F(NZ,NY,NX)
    plt_ew%CanWatP(NZ) =CanWatP(NZ,NY,NX)
    plt_pheno%WSTR(NZ) =WSTR(NZ,NY,NX)
    plt_biom%WTRVX(NZ) =WTRVX(NZ,NY,NX)
    plt_biom%CanopyLeafShethC_pft(NZ)  =CanopyLeafShethC_pft(NZ,NY,NX)

    plt_biom%ZEROL(NZ) =ZEROL(NZ,NY,NX)
    plt_biom%ZEROP(NZ) =ZEROP(NZ,NY,NX)
    plt_morph%CanopyHeight(NZ)   =CanopyHeight(NZ,NY,NX)
    plt_bgcr%ZNPP(NZ)  =ZNPP(NZ,NY,NX)


    DO L=1,NL(NY,NX)
      DO K=1,jcplx
        DO N=1,2
          DO NE=1,NumOfPlantChemElements
            plt_rbgc%RDFOME(NE,N,K,L,NZ)=RDFOME(NE,N,K,L,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DO NE=1,NumOfPlantChemElements
      DO NB=1,NumOfBranches_pft(NZ,NY,NX)
        plt_biom%EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ,NY,NX)
        plt_biom%EPOLNB(NE,NB,NZ)=EPOLNB(NE,NB,NZ,NY,NX)
        plt_biom%WTSHTBE(NE,NB,NZ)=WTSHTBE(NE,NB,NZ,NY,NX)
        plt_biom%CEPOLB(NE,NB,NZ)=CEPOLB(NE,NB,NZ,NY,NX)
        plt_biom%WTSHEBE(NE,NB,NZ)=WTSHEBE(NE,NB,NZ,NY,NX)
        plt_biom%WTSTKBE(NE,NB,NZ)=WTSTKBE(NE,NB,NZ,NY,NX)
        plt_biom%WTLFBE(NE,NB,NZ)=WTLFBE(NE,NB,NZ,NY,NX)
        plt_biom%WTRSVBE(NE,NB,NZ)=WTRSVBE(NE,NB,NZ,NY,NX)
        plt_biom%WTHSKBE(NE,NB,NZ)=WTHSKBE(NE,NB,NZ,NY,NX)
        plt_biom%WTGRBE(NE,NB,NZ)=WTGRBE(NE,NB,NZ,NY,NX)
        plt_biom%WTEARBE(NE,NB,NZ)=WTEARBE(NE,NB,NZ,NY,NX)
        plt_biom%WTNDBE(NE,NB,NZ)=WTNDBE(NE,NB,NZ,NY,NX)
        plt_biom%WGSHEXE(NE,NB,NZ)=WGSHEXE(NE,NB,NZ,NY,NX)
      ENDDO
    ENDDO

    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      plt_photo%RubiscoActivity_brpft(NB,NZ)=RubiscoActivity_brpft(NB,NZ,NY,NX)
      plt_photo%FDBKX(NB,NZ)=FDBKX(NB,NZ,NY,NX)
      plt_pheno%HourCounter4LeafOut_brch(NB,NZ)=HourCounter4LeafOut_brch(NB,NZ,NY,NX)
      plt_morph%ARLFZ(NB,NZ)=ARLFZ(NB,NZ,NY,NX)
      plt_morph%CanopyBranchLeafA_pft(NB,NZ)=CanopyBranchLeafA_pft(NB,NZ,NY,NX)

      plt_pheno%DGSTGI(NB,NZ)=DGSTGI(NB,NZ,NY,NX)
      plt_pheno%DGSTGF(NB,NZ)=DGSTGF(NB,NZ,NY,NX)
      plt_pheno%FLG4(NB,NZ)=FLG4(NB,NZ,NY,NX)
      plt_pheno%FLGZ(NB,NZ)=FLGZ(NB,NZ,NY,NX)
      plt_pheno%GROUP(NB,NZ)=GROUP(NB,NZ,NY,NX)
      plt_pheno%GSTGI(NB,NZ)=GSTGI(NB,NZ,NY,NX)
      plt_pheno%GSTGF(NB,NZ)=GSTGF(NB,NZ,NY,NX)
      plt_morph%GRNXB(NB,NZ)=GRNXB(NB,NZ,NY,NX)
      plt_morph%GRNOB(NB,NZ)=GRNOB(NB,NZ,NY,NX)
      plt_allom%GRWTB(NB,NZ)=GRWTB(NB,NZ,NY,NX)
      plt_morph%CanPBranchHeight(NB,NZ)=CanPBranchHeight(NB,NZ,NY,NX)
      plt_pheno%IDTHB(NB,NZ)=IDTHB(NB,NZ,NY,NX)
      plt_pheno%IFLGP(NB,NZ)=IFLGP(NB,NZ,NY,NX)
      plt_pheno%IFLGF(NB,NZ)=IFLGF(NB,NZ,NY,NX)
      plt_pheno%IFLGE(NB,NZ)=IFLGE(NB,NZ,NY,NX)
      plt_pheno%IFLGA(NB,NZ)=IFLGA(NB,NZ,NY,NX)
      plt_pheno%IFLGG(NB,NZ)=IFLGG(NB,NZ,NY,NX)
      plt_pheno%IFLGR(NB,NZ)=IFLGR(NB,NZ,NY,NX)
      plt_pheno%IFLGQ(NB,NZ)=IFLGQ(NB,NZ,NY,NX)
      plt_pheno%KVSTG(NB,NZ)=KVSTG(NB,NZ,NY,NX)
      plt_pheno%KVSTGN(NB,NZ)=KVSTGN(NB,NZ,NY,NX)
      plt_morph%BranchNumber_brchpft(NB,NZ)=BranchNumber_brchpft(NB,NZ,NY,NX)
      plt_morph%PSTG(NB,NZ)=PSTG(NB,NZ,NY,NX)
      plt_morph%PSTGI(NB,NZ)=PSTGI(NB,NZ,NY,NX)
      plt_morph%PSTGF(NB,NZ)=PSTGF(NB,NZ,NY,NX)
      plt_pheno%RCELX(1:NumOfPlantChemElements,NB,NZ)=RCELX(1:NumOfPlantChemElements,NB,NZ,NY,NX)
      plt_pheno%RCESX(1:NumOfPlantChemElements,NB,NZ)=RCESX(1:NumOfPlantChemElements,NB,NZ,NY,NX)
      plt_rbgc%RNH3B(NB,NZ)=RNH3B(NB,NZ,NY,NX)
      plt_pheno%TGSTGI(NB,NZ)=TGSTGI(NB,NZ,NY,NX)
      plt_pheno%TGSTGF(NB,NZ)=TGSTGF(NB,NZ,NY,NX)
      plt_pheno%VSTGX(NB,NZ)=VSTGX(NB,NZ,NY,NX)
      plt_morph%VSTG(NB,NZ)=VSTG(NB,NZ,NY,NX)
      plt_pheno%VRNY(NB,NZ)=VRNY(NB,NZ,NY,NX)
      plt_pheno%VRNZ(NB,NZ)=VRNZ(NB,NZ,NY,NX)
      plt_pheno%VRNS(NB,NZ)=VRNS(NB,NZ,NY,NX)
      plt_pheno%VRNF(NB,NZ)=VRNF(NB,NZ,NY,NX)
      plt_biom%CanPBLeafShethC(NB,NZ)=CanPBLeafShethC(NB,NZ,NY,NX)
      plt_biom%WGLFEX(1:NumOfPlantChemElements,NB,NZ) =WGLFEX(1:NumOfPlantChemElements,NB,NZ,NY,NX)
      plt_biom%WTSTXBE(1:NumOfPlantChemElements,NB,NZ)=WTSTXBE(1:NumOfPlantChemElements,NB,NZ,NY,NX)
      plt_biom%CanPBStalkC(NB,NZ)=CanPBStalkC(NB,NZ,NY,NX)
      DO M=1,pltpar%NumGrothStages
        plt_pheno%IDAY(M,NB,NZ)=IDAY(M,NB,NZ,NY,NX)
      ENDDO
      DO K=1,JNODS
        DO  L=1,JC
          DO N=1,JLI
            plt_morph%LeafA_lyrnodbrchpft(N,L,K,NB,NZ) =LeafA_lyrnodbrchpft(N,L,K,NB,NZ,NY,NX)
            plt_photo%LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ)=LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
        plt_photo%CPOOL3(K,NB,NZ)=CPOOL3(K,NB,NZ,NY,NX)
        plt_photo%CPOOL4(K,NB,NZ)=CPOOL4(K,NB,NZ,NY,NX)
        plt_photo%CO2B(K,NB,NZ)  =CO2B(K,NB,NZ,NY,NX)
        plt_photo%COMPL(K,NB,NZ) =COMPL(K,NB,NZ,NY,NX)
        plt_photo%CBXN(K,NB,NZ)  =CBXN(K,NB,NZ,NY,NX)
        plt_photo%CBXN4(K,NB,NZ) =CBXN4(K,NB,NZ,NY,NX)
        plt_photo%ETGRO(K,NB,NZ) =ETGRO(K,NB,NZ,NY,NX)
        plt_photo%ETGR4(K,NB,NZ) =ETGR4(K,NB,NZ,NY,NX)
        plt_photo%FDBK4(K,NB,NZ) =FDBK4(K,NB,NZ,NY,NX)
        plt_photo%HCOB(K,NB,NZ)  =HCOB(K,NB,NZ,NY,NX)
        plt_photo%VCGRO(K,NB,NZ) =VCGRO(K,NB,NZ,NY,NX)
        plt_photo%VGRO(K,NB,NZ)  =VGRO(K,NB,NZ,NY,NX)
        plt_photo%VCGR4(K,NB,NZ) =VCGR4(K,NB,NZ,NY,NX)
        plt_photo%VGRO4(K,NB,NZ) =VGRO4(K,NB,NZ,NY,NX)
      ENDDO
      DO K=0,JNODS
        plt_morph%ARLF1(K,NB,NZ) =ARLF(K,NB,NZ,NY,NX)
        plt_morph%HTNODX(K,NB,NZ)=HTNODX(K,NB,NZ,NY,NX)
        plt_morph%CanPSheathHeight(K,NB,NZ) =CanPSheathHeight(K,NB,NZ,NY,NX)
        plt_morph%HTNODE(K,NB,NZ)=HTNODE(K,NB,NZ,NY,NX)
        plt_biom%WGNODE(1:NumOfPlantChemElements,K,NB,NZ) =WGNODE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)
        plt_biom%WGLFE(1:NumOfPlantChemElements,K,NB,NZ)   =WGLFE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)
        plt_biom%WSLF(K,NB,NZ)   =WSLF(K,NB,NZ,NY,NX)
        plt_biom%WGSHE(1:NumOfPlantChemElements,K,NB,NZ)  =WGSHE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)
        plt_biom%WSSHE(K,NB,NZ)  =WSSHE(K,NB,NZ,NY,NX)
      ENDDO

      DO K=0,JNODS
        DO  L=1,JC
          plt_morph%CanPLNBLA(L,K,NB,NZ)=CanPLNBLA(L,K,NB,NZ,NY,NX)
          plt_biom%WGLFLE(1:NumOfPlantChemElements,L,K,NB,NZ) =WGLFLE(1:NumOfPlantChemElements,L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,JC
        plt_morph%CanopyBranchStemApft_lyr(L,NB,NZ)=CanopyBranchStemApft_lyr(L,NB,NZ,NY,NX)
      ENDDO
    enddo
    DO NE=1,NumOfPlantChemElements
      DO L=1,NL(NY,NX)
        plt_biom%WTNDLE(NE,L,NZ) =WTNDLE(NE,L,NZ,NY,NX)
      ENDDO
    ENDDO
    DO L=1,NL(NY,NX)
      DO N=1,MY(NZ,NY,NX)
        plt_biom%EPOOLR(1:NumOfPlantChemElements,N,L,NZ)=EPOOLR(1:NumOfPlantChemElements,N,L,NZ,NY,NX)
        plt_rbgc%trcs_rootml(idg_beg:idg_end-1,N,L,NZ)=trcs_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)
        plt_rbgc%trcg_rootml(idg_beg:idg_end-1,N,L,NZ)=trcg_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)

        plt_biom%RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ)=RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ,NY,NX)
        plt_biom%RootProteinConc_pftvr(N,L,NZ)=RootProteinConc_pftvr(N,L,NZ,NY,NX)

        plt_ew%PSIRoot(N,L,NZ)=PSIRoot(N,L,NZ,NY,NX)
        plt_ew%PSIRootOSMO(N,L,NZ)=PSIRootOSMO(N,L,NZ,NY,NX)
        plt_ew%PSIRootTurg(N,L,NZ)=PSIRootTurg(N,L,NZ,NY,NX)
        plt_rbgc%RUPNHB(N,L,NZ)=RUPNHB(N,L,NZ,NY,NX)
        plt_rbgc%RUPNH4(N,L,NZ)=RUPNH4(N,L,NZ,NY,NX)
        plt_rbgc%RUPH2P(N,L,NZ)=RUPH2P(N,L,NZ,NY,NX)
        plt_rbgc%RUPNOB(N,L,NZ)=RUPNOB(N,L,NZ,NY,NX)
        plt_rbgc%RUPNO3(N,L,NZ)=RUPNO3(N,L,NZ,NY,NX)
        plt_rbgc%RUPH1B(N,L,NZ)=RUPH1B(N,L,NZ,NY,NX)
        plt_rbgc%RUPH1P(N,L,NZ)=RUPH1P(N,L,NZ,NY,NX)
        plt_rbgc%RUPH2B(N,L,NZ)=RUPH2B(N,L,NZ,NY,NX)
        plt_rbgc%RUONHB(N,L,NZ)=RUONHB(N,L,NZ,NY,NX)
        plt_rbgc%RUONH4(N,L,NZ)=RUONH4(N,L,NZ,NY,NX)
        plt_rbgc%RUOH2P(N,L,NZ)=RUOH2P(N,L,NZ,NY,NX)
        plt_rbgc%RUONOB(N,L,NZ)=RUONOB(N,L,NZ,NY,NX)
        plt_rbgc%RUONO3(N,L,NZ)=RUONO3(N,L,NZ,NY,NX)
        plt_rbgc%RUOH1B(N,L,NZ)=RUOH1B(N,L,NZ,NY,NX)
        plt_rbgc%RUOH1P(N,L,NZ)=RUOH1P(N,L,NZ,NY,NX)
        plt_rbgc%RUOH2B(N,L,NZ)=RUOH2B(N,L,NZ,NY,NX)
        plt_rbgc%RUCNHB(N,L,NZ)=RUCNHB(N,L,NZ,NY,NX)
        plt_rbgc%RUCNH4(N,L,NZ)=RUCNH4(N,L,NZ,NY,NX)
        plt_rbgc%RUCH2P(N,L,NZ)=RUCH2P(N,L,NZ,NY,NX)
        plt_rbgc%RUCNOB(N,L,NZ)=RUCNOB(N,L,NZ,NY,NX)
        plt_rbgc%RUCNO3(N,L,NZ)=RUCNO3(N,L,NZ,NY,NX)
        plt_rbgc%RUCH1B(N,L,NZ)=RUCH1B(N,L,NZ,NY,NX)
        plt_rbgc%RUCH1P(N,L,NZ)=RUCH1P(N,L,NZ,NY,NX)
        plt_rbgc%RUCH2B(N,L,NZ)=RUCH2B(N,L,NZ,NY,NX)
        plt_rbgc%RCO2M(N,L,NZ)=RCO2M(N,L,NZ,NY,NX)
        plt_rbgc%RCO2N(N,L,NZ)=RCO2N(N,L,NZ,NY,NX)
        plt_rbgc%RCO2A(N,L,NZ)=RCO2A(N,L,NZ,NY,NX)
        plt_morph%PrimRootXNumL(N,L,NZ)=PrimRootXNumL(N,L,NZ,NY,NX)
        plt_morph%SecndRootXNumL(N,L,NZ)=SecndRootXNumL(N,L,NZ,NY,NX)
        plt_morph%RootLenPerP(N,L,NZ)=RootLenPerP(N,L,NZ,NY,NX)
        plt_morph%RootLenDensNLP(N,L,NZ)=RootLenDensNLP(N,L,NZ,NY,NX)
        plt_morph%RTVLP(N,L,NZ)=RTVLP(N,L,NZ,NY,NX)
        plt_morph%RTVLW(N,L,NZ)=RTVLW(N,L,NZ,NY,NX)
        plt_morph%PrimRootRadius(N,L,NZ)=PrimRootRadius(N,L,NZ,NY,NX)
        plt_morph%SecndRootRadius(N,L,NZ)=SecndRootRadius(N,L,NZ,NY,NX)
        plt_morph%RTARP(N,L,NZ)=RTARP(N,L,NZ,NY,NX)
        plt_morph%AveSecndRootLen(N,L,NZ)=AveSecndRootLen(N,L,NZ,NY,NX)
        plt_rbgc%RCO2P(N,L,NZ)=RCO2P(N,L,NZ,NY,NX)
        plt_rbgc%RUPOXP(N,L,NZ)=RUPOXP(N,L,NZ,NY,NX)
        plt_rbgc%RCO2S(N,L,NZ)=RCO2S(N,L,NZ,NY,NX)
        plt_rbgc%RUPOXS(N,L,NZ)=RUPOXS(N,L,NZ,NY,NX)
        plt_rbgc%RUPCHS(N,L,NZ)=RUPCHS(N,L,NZ,NY,NX)
        plt_rbgc%RUPN2S(N,L,NZ)=RUPN2S(N,L,NZ,NY,NX)
        plt_rbgc%RUPN3S(N,L,NZ)=RUPN3S(N,L,NZ,NY,NX)
        plt_rbgc%RUPN3B(N,L,NZ)=RUPN3B(N,L,NZ,NY,NX)
        plt_rbgc%RUPHGS(N,L,NZ)=RUPHGS(N,L,NZ,NY,NX)
        plt_rbgc%trcg_RFLA(idg_beg:idg_end-1,N,L,NZ)=trcg_RFLA(idg_beg:idg_end-1,N,L,NZ,NY,NX)
        plt_rbgc%trcg_RDFA(idg_beg:idg_end-1,N,L,NZ)=trcg_RDFA(idg_beg:idg_end-1,N,L,NZ,NY,NX)
        plt_rbgc%ROXYP(N,L,NZ) =ROXYP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNHP(N,L,NZ)=RUNNHP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNBP(N,L,NZ)=RUNNBP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNOP(N,L,NZ)=RUNNOP(N,L,NZ,NY,NX)
        plt_rbgc%RUNNXP(N,L,NZ)=RUNNXP(N,L,NZ,NY,NX)
        plt_rbgc%RUPP2P(N,L,NZ)=RUPP2P(N,L,NZ,NY,NX)
        plt_rbgc%RUPP2B(N,L,NZ)=RUPP2B(N,L,NZ,NY,NX)
        plt_rbgc%RUPP1P(N,L,NZ)=RUPP1P(N,L,NZ,NY,NX)
        plt_rbgc%RUPP1B(N,L,NZ)=RUPP1B(N,L,NZ,NY,NX)
        plt_ew%AllPlantRootH2OUptake_vr(N,L,NZ)   =AllPlantRootH2OUptake_vr(N,L,NZ,NY,NX)
        plt_rbgc%WFR(N,L,NZ)   =WFR(N,L,NZ,NY,NX)
        plt_biom%WTRTL(N,L,NZ) =WTRTL(N,L,NZ,NY,NX)
        plt_biom%PopPlantRootC_vr(N,L,NZ) =PopPlantRootC_vr(N,L,NZ,NY,NX)
        plt_biom%WSRTL(N,L,NZ) =WSRTL(N,L,NZ,NY,NX)

      enddo
      plt_bgcr%RUPNF(L,NZ) =RUPNF(L,NZ,NY,NX)
      plt_pheno%fTgrowRootP(L,NZ) =fTgrowRootP(L,NZ,NY,NX)
      plt_biom%EPOOLN(1:NumOfPlantChemElements,L,NZ)=EPOOLN(1:NumOfPlantChemElements,L,NZ,NY,NX)
    ENDDO
    DO L=1,JC
      plt_morph%CanopyStemApft_lyr(L,NZ)=CanopyStemApft_lyr(L,NZ,NY,NX)
      plt_morph%CanopyLeafApft_lyr(L,NZ)=CanopyLeafApft_lyr(L,NZ,NY,NX)
      plt_biom%CanopyLeafCpft_lyr(L,NZ) =CanopyLeafCpft_lyr(L,NZ,NY,NX)
    ENDDO
    DO N=1,MY(NZ,NY,NX)
      plt_morph%DMVL(N,NZ)=DMVL(N,NZ,NY,NX)
      plt_morph%RootPoreTortu4Gas(N,NZ)=RootPoreTortu4Gas(N,NZ,NY,NX)
      plt_morph%SecndRootXSecArea(N,NZ)=SecndRootXSecArea(N,NZ,NY,NX)
      plt_morph%PrimRootXSecArea(N,NZ)=PrimRootXSecArea(N,NZ,NY,NX)
      plt_morph%MaxPrimRootRadius1(N,NZ)=MaxPrimRootRadius1(N,NZ,NY,NX)
      plt_morph%MaxSecndRootRadius1(N,NZ)=MaxSecndRootRadius1(N,NZ,NY,NX)
      plt_morph%PrimRootSpecLen(N,NZ)=PrimRootSpecLen(N,NZ,NY,NX)
      plt_morph%RRADP(N,NZ) =RRADP(N,NZ,NY,NX)
      plt_morph%SecndRootSpecLen(N,NZ)=SecndRootSpecLen(N,NZ,NY,NX)
    ENDDO
    DO NR=1,JC
      plt_morph%NINR(NR,NZ)=NINR(NR,NZ,NY,NX)
      DO L=1,NJ(NY,NX)
        DO N=1,MY(NZ,NY,NX)
          plt_morph%PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ,NY,NX)
          plt_morph%SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ,NY,NX)
          plt_morph%RTN2(N,L,NR,NZ) =RTN2(N,L,NR,NZ,NY,NX)
          plt_biom%WTRT2E(1:NumOfPlantChemElements,N,L,NR,NZ) =WTRT2E(1:NumOfPlantChemElements,N,L,NR,NZ,NY,NX)
          plt_biom%WTRT1E(1:NumOfPlantChemElements,N,L,NR,NZ) =WTRT1E(1:NumOfPlantChemElements,N,L,NR,NZ,NY,NX)
        enddo
      enddo
      DO N=1,MY(NZ,NY,NX)
        plt_morph%PrimRootDepth(N,NR,NZ)=PrimRootDepth(N,NR,NZ,NY,NX)
        plt_biom%RTWT1E(1:NumOfPlantChemElements,N,NR,NZ)=RTWT1E(1:NumOfPlantChemElements,N,NR,NZ,NY,NX)
      enddo
    enddo
    DO NE=1,NumOfPlantChemElements
      DO M=1,jsken
        plt_biom%WTSTDE(NE,M,NZ)=WTSTDE(NE,M,NZ,NY,NX)
      ENDDO
    ENDDO
    DO L=0,NJ(NY,NX)
      DO K=1,micpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          plt_bgcr%ESNC(1:NumOfPlantChemElements,M,K,L,NZ)=ESNC(1:NumOfPlantChemElements,M,K,L,NZ,NY,NX)
        enddo
      enddo
    ENDDO
    DO NE=1,NumOfPlantChemElements
      DO M=1,jsken
        DO N=0,JP
          plt_soilchem%CFOPE(NE,N,M,NZ)=CFOPE(NE,N,M,NZ,NY,NX)
        enddo
      enddo
    ENDDO
  ENDDO

  DO L=1,NL(NY,NX)
    DO M=1,NPH
      plt_rbgc%ROXSK(M,L)=ROXSK(M,L,NY,NX)
      plt_soilchem%THETPM(M,L)=THETPM(M,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine PlantAPISend


!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMSend(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,N,M,NN,NZ,K,NB
!  Integers
  plt_site%IETYP=IETYP(NY,NX)
  plt_site%NP=NP(NY,NX)
  plt_site%NU=NU(NY,NX)
  plt_morph%StemAreag=StemAreag(NY,NX)
  plt_morph%CanopyLA_grd=CanopyLA_grd(NY,NX)
  plt_site%ZEROS =ZEROS(NY,NX)
  plt_site%ZERO  =ZERO
  plt_ew%SnowDepth   =SnowDepth(NY,NX)
  plt_ew%TairK     =TairK(NY,NX)
  plt_morph%GridMaxCanopyHeight   =GridMaxCanopyHeight(NY,NX)
  plt_site%WindMesHeight=WindMesHeight(NY,NX)
  plt_ew%ZeroPlanDisp=ZeroPlanDisp(NY,NX)
  plt_site%WindSpeedAtm=WindSpeedAtm(NY,NX)
  plt_ew%VLHeatCapSnowMin=VLHeatCapSnowMin(NY,NX)
  plt_ew%VLHeatCapSurfSnow=VLHeatCapSnow(1,NY,NX)
  plt_morph%CanopyHeightz(0)=CanopyHeightz(0,NY,NX)
  DO L=1,JC
    plt_morph%CanopyStemA_lyr(L)=CanopyStemA_lyr(L,NY,NX)
    plt_morph%CanopyLAgrid_lyr(L)=CanopyLAgrid_lyr(L,NY,NX)
    plt_morph%CanopyHeightz(L)=CanopyHeightz(L,NY,NX)
    plt_rad%TAUS(L)=TAUS(L,NY,NX)
  ENDDO
  plt_rad%TAUS(JC+1)=TAUS(JC+1,NY,NX)

  DO L=0,NL(NY,NX)
    plt_site%AREA3(L)=AREA(3,L,NY,NX)
    plt_soilchem%VLSoilPoreMicP(L)=VLSoilPoreMicP(L,NY,NX)
    plt_soilchem%VLSoilMicP(L)=VLSoilMicP(L,NY,NX)
    plt_soilchem%VLWatMicP(L)=VLWatMicP(L,NY,NX)
  ENDDO
  plt_rad%SSIN=SSIN(NY,NX)
  plt_site%ZNOON=ZNOON(NY,NX)
  plt_morph%CanopyArea_grid=CanopyArea_grid(NY,NX)
  plt_rad%GAZI=GAZI(NY,NX)
  plt_rad%GCOS=GCOS(NY,NX)
  plt_rad%GSIN=GSIN(NY,NX)
  plt_rad%RADY=RADY(NY,NX)
  plt_rad%RAPY=RAPY(NY,NX)
  plt_rad%RADS=RADS(NY,NX)
  plt_rad%RAPS=RAPS(NY,NX)
  plt_site%SoiSurfRoughnesst0=SoiSurfRoughnesst0(NY,NX)
  plt_ew%VcumWatSnow=VcumWatSnow(NY,NX)
  plt_ew%VcumIceSnow=VcumIceSnow(NY,NX)
  plt_ew%VcumDrySnoWE=VcumDrySnoWE(NY,NX)
  plt_rad%TYSIN  =TYSIN
  plt_rad%SoilAlbedo   =SoilAlbedo(NY,NX)
  plt_rad%ALBX   =ALBX(NY,NX)
  plt_site%ZEROS2=ZEROS2(NY,NX)
  plt_site%POROS1=POROS(NU(NY,NX),NY,NX)
  DO NZ=1,NP(NY,NX)
    plt_morph%CanopyLeafA_pft(NZ)=CanopyLeafA_pft(NZ,NY,NX)
    plt_morph%CanopyHeight(NZ)=CanopyHeight(NZ,NY,NX)
    plt_morph%ClumpFactort(NZ)=ClumpFactort(NZ,NY,NX)
    plt_rad%CanopySWabsorpty_pft(NZ)=CanopySWabsorpty_pft(NZ,NY,NX)
    plt_rad%CanopyPARabsorpty_pft(NZ)=CanopyPARabsorpty_pft(NZ,NY,NX)
    plt_rad%TAUR(NZ)=TAUR(NZ,NY,NX)
    plt_rad%CanopySWAlbedo_pft(NZ)=CanopySWAlbedo_pft(NZ,NY,NX)
    plt_rad%TAUP(NZ)=TAUP(NZ,NY,NX)
    plt_rad%CanopyPARalbedo_pft(NZ)=CanopyPARalbedo_pft(NZ,NY,NX)
    plt_morph%NumOfBranches_pft(NZ)=NumOfBranches_pft(NZ,NY,NX)
    plt_morph%ClumpFactor(NZ)=ClumpFactor(NZ,NY,NX)
    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO K=0,JNODS
        DO  L=1,JC
          plt_morph%CanPLNBLA(L,K,NB,NZ)=CanPLNBLA(L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,JC
        plt_morph%CanopyBranchStemApft_lyr(L,NB,NZ)=CanopyBranchStemApft_lyr(L,NB,NZ,NY,NX)
      ENDDO
      DO K=1,JNODS
        DO  L=1,JC
          DO N=1,JLI
            plt_morph%LeafA_lyrnodbrchpft(N,L,K,NB,NZ)=LeafA_lyrnodbrchpft(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
      DO  L=1,JC
        DO N=1,JLI
          plt_morph%StemA_lyrnodbrchpft(N,L,NB,NZ)=StemA_lyrnodbrchpft(N,L,NB,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO N=1,JSA
    plt_rad%OMEGAG(N)=OMEGAG(N,NY,NX)
  ENDDO
  DO N=1,JLI
    plt_rad%ZCOS(N)=ZCOS(N)
    plt_rad%ZSIN(N)=ZSIN(N)
  ENDDO
  DO NN=1,JLA
    DO M=1,JLI
      DO N=1,JSA
        plt_rad%OMEGA(N,M,NN)=OMEGA(N,M,NN)
        plt_rad%OMEGX(N,M,NN)=OMEGX(N,M,NN)
        plt_rad%IALBY(N,M,NN)=IALBY(N,M,NN)
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMSend

!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMRecv(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  integer :: N,M,NN,L,NZ,K,NB

  ZeroPlanDisp(NY,NX)=plt_ew%ZeroPlanDisp
  RoughHeight(NY,NX)=plt_ew%RoughHeight
  BndlResistAboveCanG(NY,NX)=plt_ew%BndlResistAboveCanG
  RIB(NY,NX)=plt_ew%RIB
  GridMaxCanopyHeight(NY,NX)=plt_morph%GridMaxCanopyHeight
  RADS(NY,NX)=plt_rad%RADS
  RADY(NY,NX)=plt_rad%RADY
  RAPS(NY,NX)=plt_rad%RAPS
  RAPY(NY,NX)=plt_rad%RAPY
  SWRadOnGrnd(NY,NX)=plt_rad%SWRadOnGrnd
  FracSWRad2Grnd(NY,NX)=plt_rad%FracSWRad2Grnd
  RAD(NY,NX)=plt_rad%RAD0
  RAP(NY,NX)=plt_rad%RAP0
  DO L=0,JC
    CanopyHeightz(L,NY,NX)=plt_morph%CanopyHeightz(L)
  ENDDO
  DO L=1,JC
    TAUS(L,NY,NX)=plt_rad%TAUS(L)
    TAU0(L,NY,NX)=plt_rad%TAU0(L)
  ENDDO
  CanopyArea_grid(NY,NX)=plt_morph%CanopyArea_grid
  DO NZ=1,NP(NY,NX)
    CanopyArea_pft(NZ,NY,NX)=plt_morph%CanopyArea_pft(NZ)
    SWRadByCanP(NZ,NY,NX) =plt_rad%SWRadByCanP(NZ)
    PARByCanP(NZ,NY,NX) =plt_rad%PARByCanP(NZ)
    ClumpFactort(NZ,NY,NX)  =plt_morph%ClumpFactort(NZ)
    FracPARByCanP(NZ,NY,NX)=plt_rad%FracPARByCanP(NZ)

    DO L=1,JC
      DO M=1,JSA
        DO  N=1,JLI
          PARDIF(N,M,L,NZ,NY,NX)=plt_rad%PARDIF(N,M,L,NZ)
          PAR(N,M,L,NZ,NY,NX)   =plt_rad%PAR(N,M,L,NZ)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMRecv

end module PlantAPI
