module NoduleBGCMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : safe_adb,AZMAX1
  use EcosimConst
  use PlantAPIData
  use PlantMathFuncMod
  use GrosubPars
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: CanopyNoduleBiochemistry
  public :: RootNodulBiochemistry
  contains

!------------------------------------------------------------------------------------------

  subroutine CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WFNG,CanopyN2Fix_pft)
  !
  !Nodule update for branch NB on pft NZ
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: TFN5,WFNG
  real(r8), intent(inout) :: CanopyN2Fix_pft(JP1)
  integer :: M
  real(r8) :: ZPOOLD
  real(r8) :: CCC,CNC,CPC
  REAL(R8) :: cpoolt
  real(r8) :: CPOOLD
  real(r8) :: NodulNonstElmConc(NumPlantChemElms)
  real(r8) :: NodulUseNonstElm_brch(NumPlantChemElms)
  real(r8) :: CCNDLB
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: NoduleBiomCGrowth
  real(r8) :: PPOOLD
  real(r8) :: RCO2T
  real(r8) :: RespNonst_Oltd,Rmaint,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,Rauto4Nfix_brch
  real(r8) :: RCanopyN2Fix
  real(r8) :: NoduleElmDecayLoss(NumPlantChemElms)
  real(r8) :: NoduleElmntDecay2Litr(NumPlantChemElms)
  real(r8) :: NodulElmDecayRecyc(NumPlantChemElms)
  real(r8) :: NoduleCResp
  real(r8) :: NodulELmLoss2Senes(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Litr(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Recyc(NumPlantChemElms)
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTLSB1,WTNDB1,WTLSBT
  real(r8) :: XFRE(NumPlantChemElms)
  REAL(R8) :: RCCC,RCCN,RCCP
  integer :: NE
!     begin_execution
  associate(                        &
    NU                       =>  plt_site%NU       , &
    ZERO                     =>  plt_site%ZERO     , &
    AREA3                    =>  plt_site%AREA3    , &
    k_fine_litr              => pltpar%k_fine_litr,&
    CFOPE                    =>  plt_soilchem%CFOPE, &
    iPlantNfixType           =>  plt_morph%iPlantNfixType  , &
    fTgrowCanP               =>  plt_pheno%fTgrowCanP    , &
    GrossResp_pft            =>  plt_bgcr%GrossResp_pft    , &
    ECO_ER_col               =>  plt_bgcr%ECO_ER_col     , &
    CanopyPlusNoduRespC_pft  =>  plt_bgcr%CanopyPlusNoduRespC_pft    , &
    Eco_AutoR_col            =>  plt_bgcr%Eco_AutoR_col     , &
    CO2NetFix_pft            =>  plt_bgcr%CO2NetFix_pft     , &
    LitfalChemElm_pvr        =>  plt_bgcr%LitfalChemElm_pvr     , &
    ifoliar                  =>  pltpar%ifoliar    , &
    NoduGrowthYield_pft      =>  plt_allom%NoduGrowthYield_pft   , &
    NodulerNC_pft            =>  plt_allom%NodulerNC_pft   , &
    NodulerPC_pft            =>  plt_allom%NodulerPC_pft   , &
    LeafPetolBiomassC_brch   =>  plt_biom%LeafPetolBiomassC_brch    , &
    NonstructElm_brch        =>  plt_biom%NonstructElm_brch   , &
    NodulNonstElm_brch       =>  plt_biom%NodulNonstElm_brch   , &
    ZEROP                    =>  plt_biom%ZEROP    , &
    ZEROL                    =>  plt_biom%ZEROL    , &
    CanopyNodulChemElm_brch =>  plt_biom%CanopyNodulChemElm_brch     &
  )
!     iPlantNfixType=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
!
  IF(is_canopy_N2fix(iPlantNfixType(NZ)))THEN
!
!     INITIAL INFECTION
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     NodulBiomCatInfection=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,NodulerPC_pft=bacterial N:C,P:C ratio from PFT file
!
    IF(CanopyNodulChemElm_brch(ielmc,NB,NZ).LE.0.0_r8)THEN
      CanopyNodulChemElm_brch(ielmc,NB,NZ)=CanopyNodulChemElm_brch(ielmc,NB,NZ)+NodulBiomCatInfection*AREA3(NU)
      CanopyNodulChemElm_brch(ielmn,NB,NZ)=CanopyNodulChemElm_brch(ielmn,NB,NZ)+NodulBiomCatInfection*AREA3(NU)*NodulerNC_pft(NZ)
      CanopyNodulChemElm_brch(ielmp,NB,NZ)=CanopyNodulChemElm_brch(ielmp,NB,NZ)+NodulBiomCatInfection*AREA3(NU)*NodulerPC_pft(NZ)
    ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
    IF(CanopyNodulChemElm_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      NodulNonstElmConc(ielmc)=AZMAX1(NodulNonstElm_brch(ielmc,NB,NZ)/CanopyNodulChemElm_brch(ielmc,NB,NZ))
      NodulNonstElmConc(ielmn)=AZMAX1(NodulNonstElm_brch(ielmn,NB,NZ)/CanopyNodulChemElm_brch(ielmc,NB,NZ))
      NodulNonstElmConc(ielmp)=AZMAX1(NodulNonstElm_brch(ielmp,NB,NZ)/CanopyNodulChemElm_brch(ielmc,NB,NZ))
    ELSE
      NodulNonstElmConc(ielmc)=1.0_r8
      NodulNonstElmConc(ielmn)=1.0_r8
      NodulNonstElmConc(ielmp)=1.0_r8
    ENDIF

    IF(NodulNonstElmConc(ielmc).GT.ZERO)THEN
      CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmn)+NodulNonstElmConc(ielmc)*CNKI) &
        ,safe_adb(NodulNonstElmConc(ielmp),NodulNonstElmConc(ielmp)+NodulNonstElmConc(ielmc)*CPKI)))
      CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmn)/CNKI)))
      CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmp)/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    IF(CanopyNodulChemElm_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      FCNPF=AMIN1(1.0_r8 &
        ,SQRT(CanopyNodulChemElm_brch(ielmn,NB,NZ)/(CanopyNodulChemElm_brch(ielmc,NB,NZ)*NodulerNC_pft(NZ))) &
        ,SQRT(CanopyNodulChemElm_brch(ielmp,NB,NZ)/(CanopyNodulChemElm_brch(ielmc,NB,NZ)*NodulerPC_pft(NZ))))
    ELSE
      FCNPF=1.0_r8
    ENDIF
    !MM factor for nodule maintenance respiration
    SPNDLI=NodulNonstElmConc(ielmc)/(NodulNonstElmConc(ielmc)+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RespNonst_Oltd=respiration from non-structural C
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDB=bacterial C mass
!     fTgrowCanP=temperature function for canopy growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNG=growth function of canopy water potential
!
    RespNonst_Oltd=AZMAX1(AMIN1(NodulNonstElm_brch(ielmc,NB,NZ),&
      VMXO*CanopyNodulChemElm_brch(ielmc,NB,NZ))*FCNPF*fTgrowCanP(NZ)*WFNG)
!     CPOOLNX=NodulNonstElm_brch(ielmc,NB,NZ)
!     VMXOX=VMXO*CanopyNodulChemElm_brch(ielmc,NB,NZ)*FCNPF*fTgrowCanP(NZ)*WFNG
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     Rmaint=bacterial maintenance respiration
!     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN5=temperature function for canopy maintenance respiration
!     WTNDBN=bacterial N mass
!
    Rmaint=AZMAX1(RmSpecPlant*TFN5*CanopyNodulChemElm_brch(ielmn,NB,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDL=difference between non-structural C respn and mntc respn
!     RGNDL=growth respiration unlimited by N,P
!     RSNDL=excess maintenance respiration
!
    RXNDL=RespNonst_Oltd-Rmaint
    RGNDL=AZMAX1(RXNDL)
    RSNDL=AZMAX1(-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDB,WTNDBN=bacterial C,N mass
!     NodulerNC_pft=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     Rauto4Nfix_brch=respiration for N2 fixation
!     RCanopyN2Fix,CanopyN2Fix_pft=branch,total N2 fixation
!
    RGN2P=AZMAX1(CanopyNodulChemElm_brch(ielmc,NB,NZ)*NodulerNC_pft(NZ)-&
      CanopyNodulChemElm_brch(ielmn,NB,NZ))/EN2F

    IF(RGNDL.GT.ZEROP(NZ))THEN
      Rauto4Nfix_brch=RGNDL*RGN2P/(RGNDL+RGN2P)
    ELSE
      Rauto4Nfix_brch=0._r8
    ENDIF
    RCanopyN2Fix=Rauto4Nfix_brch*EN2F
    CanopyN2Fix_pft(NZ)=CanopyN2Fix_pft(NZ)+RCanopyN2Fix
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION AND LEAKAGE
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     LeafPetolBiomassC_brch=leaf+petiole mass
!     CCNDLB=bacteria:leaf+petiole ratio
!     RDNDBX=effect of CCNDLB on bacteria decomposition rate
!     SPNDX=specific bacterial decomposition rate at current CCNDLB
!     SPNDL=specific decomposition rate by bacterial N2 fixers
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     NoduleElmDecayLoss(ielmc),NoduleElmDecayLoss(ielmn),NoduleElmDecayLoss(ielmp)=bacterial C,N,P loss from decomposition
!     NoduleElmntDecay2Litr(ielmc),NoduleElmntDecay2Litr(ielmn),NoduleElmntDecay2Litr(ielmp)=bacterial C,N,P decomposition to LitrFall
!     NodulElmDecayRecyc(ielmc),NodulElmDecayRecyc(ielmn),NodulElmDecayRecyc(ielmp)=bacterial C,N,P decomposition to recycling
!
    RCCC=RCCZN+CCC*RCCYN
    RCCN=CNC*RCCXN
    RCCP=CPC*RCCQN
    SPNDX=SPNDL*SQRT(fTgrowCanP(NZ)*WFNG)
    DO NE=1,NumPlantChemElms
      NoduleElmDecayLoss(NE)=SPNDX*CanopyNodulChemElm_brch(NE,NB,NZ)
    ENDDO

    NoduleElmntDecay2Litr(ielmc)=NoduleElmDecayLoss(ielmc)*(1.0_r8-RCCC)
    NoduleElmntDecay2Litr(ielmn)=NoduleElmDecayLoss(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    NoduleElmntDecay2Litr(ielmp)=NoduleElmDecayLoss(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)

    DO NE=1,NumPlantChemElms
      NodulElmDecayRecyc(NE)=NoduleElmDecayLoss(NE)-NoduleElmntDecay2Litr(NE)
    ENDDO
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     NodulUseNonstElm_brch(ielmc)=total non-structural C used in bacterial growth and growth respiration
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     Rmaint=bacterial maintenance respiration
!     RespNonst_Oltd=respiration from non-structural C
!     NodulElmDecayRecyc(ielmc)=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     Rauto4Nfix_brch=respiration for N2 fixation
!     NoduleBiomCGrowth=bacterial growth
!     NoduGrowthYield_pft=bacterial growth yield
!     NoduleCResp=bacterial respiration for growth and N2 fixation
!     NodulUseNonstElm_brch(ielmn),NodulUseNonstElm_brch(ielmp)=nonstructural N,P used in growth
!     CNND,NodulerPC_pft=bacterial N:C,P:C ratio from PFT file
!     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
    NodulUseNonstElm_brch(ielmc)=AMIN1(NodulNonstElm_brch(ielmc,NB,NZ)-AMIN1(Rmaint,RespNonst_Oltd) &
      -Rauto4Nfix_brch+NodulElmDecayRecyc(ielmc),(RGNDL-Rauto4Nfix_brch)/(1.0_r8-NoduGrowthYield_pft(NZ)))
    NoduleBiomCGrowth=NodulUseNonstElm_brch(ielmc)*NoduGrowthYield_pft(NZ)
    NoduleCResp=Rauto4Nfix_brch+NodulUseNonstElm_brch(ielmc)*(1.0_r8-NoduGrowthYield_pft(NZ))
    NodulUseNonstElm_brch(ielmn)=AZMAX1(AMIN1(NodulNonstElm_brch(ielmn,NB,NZ),NoduleBiomCGrowth*NodulerNC_pft(NZ))) &
      *NodulNonstElmConc(ielmn)/(NodulNonstElmConc(ielmn)+CZKM)
    NodulUseNonstElm_brch(ielmp)=AZMAX1(AMIN1(NodulNonstElm_brch(ielmp,NB,NZ),NoduleBiomCGrowth*NodulerPC_pft(NZ))) &
      *NodulNonstElmConc(ielmp)/(NodulNonstElmConc(ielmp)+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmp)=bacterial C,N,P loss from senescence
!     NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmp)=bacterial C,N,P senescence to LitrFall
!     NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmp)=bacterial C,N,P senescence to recycling
!
    IF(RSNDL.GT.0.0_r8 .AND. CanopyNodulChemElm_brch(ielmc,NB,NZ).GT.ZEROP(NZ) .AND. RCCC.GT.ZERO)THEN
      NodulELmLoss2Senes(ielmc)=RSNDL/RCCC
      NodulELmLoss2Senes(ielmn)=NodulELmLoss2Senes(ielmc)*CanopyNodulChemElm_brch(ielmn,NB,NZ)/CanopyNodulChemElm_brch(ielmc,NB,NZ)
      NodulELmLoss2Senes(ielmp)=NodulELmLoss2Senes(ielmc)*CanopyNodulChemElm_brch(ielmp,NB,NZ)/CanopyNodulChemElm_brch(ielmc,NB,NZ)

      NodulELmSenes2Litr(ielmc)=NodulELmLoss2Senes(ielmc)*(1.0_r8-RCCC)
      NodulELmSenes2Litr(ielmn)=NodulELmLoss2Senes(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
      NodulELmSenes2Litr(ielmp)=NodulELmLoss2Senes(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
      DO NE=1,NumPlantChemElms
        NodulELmSenes2Recyc(NE)=NodulELmLoss2Senes(NE)-NodulELmSenes2Litr(NE)
      ENDDO
    ELSE
      NodulELmLoss2Senes(1:NumPlantChemElms)=0._r8
      NodulELmSenes2Litr(1:NumPlantChemElms)=0._r8
      NodulELmSenes2Recyc(1:NumPlantChemElms)=0._r8
    ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2T=total C respiration
!     Rmaint=bacterial maintenance respiration
!     RespNonst_Oltd=respiration from non-structural C
!     NoduleCResp=bacterial respiration for growth and N2 fixation
!     NodulELmSenes2Recyc(ielmc)=bacterial C senescence to recycling
!     GrossResp_pft,CanopyPlusNoduRespC_pft=total,above-ground PFT respiration
!     CO2NetFix_pft=PFT net CO2 fixation
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!
    RCO2T=AMIN1(Rmaint,RespNonst_Oltd)+NoduleCResp+NodulELmSenes2Recyc(ielmc)
    GrossResp_pft(NZ)=GrossResp_pft(NZ)-RCO2T
    CanopyPlusNoduRespC_pft(NZ)=CanopyPlusNoduRespC_pft(NZ)-RCO2T
    CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-RCO2T
    ECO_ER_col=ECO_ER_col-RCO2T
    Eco_AutoR_col=Eco_AutoR_col-RCO2T
!
!     NODULE LitrFall CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     NoduleElmntDecay2Litr(ielmc),NoduleElmntDecay2Litr(ielmn),NoduleElmntDecay2Litr(ielmp)=bacterial C,N,P decomposition to LitrFall
!     NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmp)=bacterial C,N,P senescence to LitrFall
!
    D6470: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        LitfalChemElm_pvr(NE,M,k_fine_litr,0,NZ)=LitfalChemElm_pvr(NE,M,k_fine_litr,0,NZ) &
          +CFOPE(NE,ifoliar,M,NZ)*(NoduleElmntDecay2Litr(NE)+NodulELmSenes2Litr(NE))
      ENDDO
    ENDDO D6470
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     Rmaint=bacterial maintenance respiration
!     RespNonst_Oltd=respiration from non-structural C
!     Rauto4Nfix_brch=respiration for N2 fixation
!     NodulUseNonstElm_brch(ielmc)=total non-structural C used in bacterial growth and growth respiration
!     NodulElmDecayRecyc(ielmc),NodulElmDecayRecyc(ielmn),NodulElmDecayRecyc(ielmp)=bacterial C,N,P decomposition to recycling
!     NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmp)=bacterial C,N,P senescence to recycling
!     NodulUseNonstElm_brch(ielmn),NodulUseNonstElm_brch(ielmp)=nonstructural N,P used in growth
!     RCanopyN2Fix=branch N2 fixation
!
    NodulNonstElm_brch(ielmc,NB,NZ)=NodulNonstElm_brch(ielmc,NB,NZ)-AMIN1(Rmaint,RespNonst_Oltd)-Rauto4Nfix_brch&
      -NodulUseNonstElm_brch(ielmc)+NodulElmDecayRecyc(ielmc)
    NodulNonstElm_brch(ielmn,NB,NZ)=NodulNonstElm_brch(ielmn,NB,NZ) &
      -NodulUseNonstElm_brch(ielmn)+NodulElmDecayRecyc(ielmn)+NodulELmSenes2Recyc(ielmn)+RCanopyN2Fix
    NodulNonstElm_brch(ielmp,NB,NZ)=NodulNonstElm_brch(ielmp,NB,NZ) &
      -NodulUseNonstElm_brch(ielmp)+NodulElmDecayRecyc(ielmp)+NodulELmSenes2Recyc(ielmp)
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     NoduleBiomCGrowth=bacterial growth
!     NoduleElmDecayLoss(ielmc),NoduleElmDecayLoss(ielmn),NoduleElmDecayLoss(ielmp)=bacterial C,N,P loss from decomposition
!     NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmp)=bacterial C,N,P loss from senescence
!     NodulUseNonstElm_brch(ielmn),NodulUseNonstElm_brch(ielmp)=nonstructural N,P used in growth
!
    CanopyNodulChemElm_brch(ielmc,NB,NZ)=CanopyNodulChemElm_brch(ielmc,NB,NZ)+NoduleBiomCGrowth &
      -NoduleElmDecayLoss(ielmc)-NodulELmLoss2Senes(ielmc)
    CanopyNodulChemElm_brch(ielmn,NB,NZ)=CanopyNodulChemElm_brch(ielmn,NB,NZ) &
      +NodulUseNonstElm_brch(ielmn)-NoduleElmDecayLoss(ielmn)-NodulELmLoss2Senes(ielmn)
    CanopyNodulChemElm_brch(ielmp,NB,NZ)=CanopyNodulChemElm_brch(ielmp,NB,NZ) &
      +NodulUseNonstElm_brch(ielmp)-NoduleElmDecayLoss(ielmp)-NodulELmLoss2Senes(ielmp)
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN BRANCH AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     LeafPetolBiomassC_brch=leaf+petiole C mass
!     WTNDB=bacterial C mass
!     NodulBiomCatInfection=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGB=parameter to calculate nonstructural C,N,P exchange
!     CCNDLB=bacteria:leaf+petiole ratio
!     XFRE(ielmc),XFRE(ielmn),XFRE(ielmp)=nonstructural C,N,P transfer
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!
    IF(NonstructElm_brch(ielmc,NB,NZ).GT.ZEROP(NZ).AND.LeafPetolBiomassC_brch(NB,NZ).GT.ZEROL(NZ))THEN
      CCNDLB=CanopyNodulChemElm_brch(ielmc,NB,NZ)/LeafPetolBiomassC_brch(NB,NZ)
      WTLSB1=LeafPetolBiomassC_brch(NB,NZ)
      WTNDB1=AMIN1(LeafPetolBiomassC_brch(NB,NZ),AMAX1(NodulBiomCatInfection*AREA3(NU) &
        ,CanopyNodulChemElm_brch(ielmc,NB,NZ)))
      WTLSBT=WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZEROP(NZ))THEN
        FXRNX=FXRN(iPlantNfixType(NZ))/(1.0_r8+CCNDLB/CCNGB)
!    2/(1.0+CCNDLB/(CCNGB*FXRN(iPlantNfixType(NZ))))
        CPOOLD=(NonstructElm_brch(ielmc,NB,NZ)*WTNDB1-NodulNonstElm_brch(ielmc,NB,NZ)*WTLSB1)/WTLSBT
        XFRE(ielmc)=FXRNX*CPOOLD
        NonstructElm_brch(ielmc,NB,NZ)=NonstructElm_brch(ielmc,NB,NZ)-XFRE(ielmc)
        NodulNonstElm_brch(ielmc,NB,NZ)=NodulNonstElm_brch(ielmc,NB,NZ)+XFRE(ielmc)

        CPOOLT=NonstructElm_brch(ielmc,NB,NZ)+NodulNonstElm_brch(ielmc,NB,NZ)
        IF(CPOOLT.GT.ZEROP(NZ))THEN
          ZPOOLD=(NonstructElm_brch(ielmn,NB,NZ)*NodulNonstElm_brch(ielmc,NB,NZ)-NodulNonstElm_brch(ielmn,NB,NZ)*NonstructElm_brch(ielmc,NB,NZ))/CPOOLT
          XFRE(ielmn)=FXRNX*ZPOOLD
          PPOOLD=(NonstructElm_brch(ielmp,NB,NZ)*NodulNonstElm_brch(ielmc,NB,NZ)-NodulNonstElm_brch(ielmp,NB,NZ)*NonstructElm_brch(ielmc,NB,NZ))/CPOOLT
          XFRE(ielmp)=FXRNX*PPOOLD
          DO NE=2,NumPlantChemElms
            NonstructElm_brch(NE,NB,NZ)=NonstructElm_brch(NE,NB,NZ)-XFRE(NE)
            NodulNonstElm_brch(NE,NB,NZ)=NodulNonstElm_brch(NE,NB,NZ)+XFRE(NE)
          ENDDO
        ENDIF

      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine CanopyNoduleBiochemistry

!------------------------------------------------------------------------------------------

  subroutine RootNodulBiochemistry(I,J,NZ,TFN6_vr,fRootGrowPSISense)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: fRootGrowPSISense(2,JZ1)
  integer :: L,M
  real(r8) :: ZPOOLD
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLT
  real(r8) :: CPOOLD
  real(r8) :: NodulNonstElmConc(NumPlantChemElms)
  real(r8) :: NodulUseNonstElm_rvr(NumPlantChemElms),CPOOLNX
  real(r8) :: CCNDLR
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: NoduleBiomCGrowth
  real(r8) :: PPOOLD
  real(r8) :: RCO2T
  real(r8) :: RCO2TM
  real(r8) :: RespNonst_Oltd,Rmaint,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,Rauto4Nfix_rvr
  real(r8) :: NoduleElmDecayLoss(NumPlantChemElms)
  real(r8) :: NoduleElmntDecay2Litr(NumPlantChemElms)
  real(r8) :: NodulElmDecayRecyc(NumPlantChemElms)
  real(r8) :: NoduleCResp
  real(r8) :: NodulELmLoss2Senes(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Litr(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Recyc(NumPlantChemElms)
  real(r8) :: RespNonst_OUltd,RXNDLM,RGNDLM
  real(r8) :: RSNDLM
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTRTD1,WTNDL1,WTRTDT
  real(r8) :: XFRE(NumPlantChemElms)
  real(r8) :: RCCC,RCCN,RCCP
  integer  :: NE
!     begin_execution
  associate(                          &
    NU                           =>   plt_site%NU        , &
    AREA3                        =>   plt_site%AREA3     , &
    ZERO                         =>   plt_site%ZERO      , &
    fTgrowRootP_vr               =>   plt_pheno%fTgrowRootP_vr     , &
    NoduGrowthYield_pft          =>   plt_allom%NoduGrowthYield_pft    , &
    NodulerNC_pft                =>   plt_allom%NodulerNC_pft    , &
    NodulerPC_pft                =>   plt_allom%NodulerPC_pft    , &
    k_fine_litr                  =>   pltpar%k_fine_litr , &
    iroot                        =>   pltpar%iroot       , &
    RootRespPotent_pvr           =>   plt_rbgc%RootRespPotent_pvr     , &
    RCO2N_pvr                    =>   plt_rbgc%RCO2N_pvr     , &
    RootAutoRO2Limiter_pvr       =>   plt_rbgc%RootAutoRO2Limiter_pvr      , &
    RCO2A_pvr                    =>   plt_rbgc%RCO2A_pvr     , &
    LitfalChemElm_pvr            =>   plt_bgcr%LitfalChemElm_pvr      , &
    RootN2Fix_pft                =>   plt_rbgc%RootN2Fix_pft      , &
    RootN2Fix_pvr                =>   plt_bgcr%RootN2Fix_pvr     , &
    PopuPlantRootC_vr            =>   plt_biom% PopuPlantRootC_vr    , &
    RootNodulElm_pvr             =>   plt_biom%RootNodulElm_pvr    , &
    ZEROP                        =>   plt_biom%ZEROP     , &
    RootNodulNonstElm_pvr        =>   plt_biom%RootNodulNonstElm_pvr   , &
    ZEROL                        =>   plt_biom%ZEROL     , &
    RootMycoNonstElm_pvr         =>   plt_biom%RootMycoNonstElm_pvr    , &
    CFOPE                        =>   plt_soilchem%CFOPE , &
    iPlantNfixType               =>   plt_morph%iPlantNfixType   , &
    NIXBotRootLayer_pft          =>   plt_morph%NIXBotRootLayer_pft       &
  )
!     iPlantNfixType=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     NodulBiomCatInfection=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,NodulerPC_pft=bacterial N:C,P:C ratio from PFT file
!
  IF(is_root_N2fix(iPlantNfixType(NZ)))THEN
    D5400: DO L=NU,NIXBotRootLayer_pft(NZ)
      IF(PopuPlantRootC_vr(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
!
!     INITIAL INFECTION
!
        IF(RootNodulElm_pvr(ielmc,L,NZ).LE.0.0_r8)THEN
          RootNodulElm_pvr(ielmc,L,NZ)=RootNodulElm_pvr(ielmc,L,NZ)+NodulBiomCatInfection*AREA3(NU)
          RootNodulElm_pvr(ielmn,L,NZ)=RootNodulElm_pvr(ielmn,L,NZ)+NodulBiomCatInfection*AREA3(NU)*NodulerNC_pft(NZ)
          RootNodulElm_pvr(ielmp,L,NZ)=RootNodulElm_pvr(ielmp,L,NZ)+NodulBiomCatInfection*AREA3(NU)*NodulerPC_pft(NZ)
        ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
        IF(RootNodulElm_pvr(ielmc,L,NZ).GT.ZEROP(NZ))THEN
          NodulNonstElmConc(ielmc)=AZMAX1(RootNodulNonstElm_pvr(ielmc,L,NZ)/RootNodulElm_pvr(ielmc,L,NZ))
          NodulNonstElmConc(ielmn)=AZMAX1(RootNodulNonstElm_pvr(ielmn,L,NZ)/RootNodulElm_pvr(ielmc,L,NZ))
          NodulNonstElmConc(ielmp)=AZMAX1(RootNodulNonstElm_pvr(ielmp,L,NZ)/RootNodulElm_pvr(ielmc,L,NZ))
        ELSE
          NodulNonstElmConc(ielmc)=1.0_r8
          NodulNonstElmConc(ielmn)=1.0_r8
          NodulNonstElmConc(ielmp)=1.0_r8
        ENDIF
        IF(NodulNonstElmConc(ielmc).GT.ZERO)THEN
          CCC=AZMAX1(AMIN1(1.0_r8 &
            ,safe_adb(NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmn)+NodulNonstElmConc(ielmc)*CNKI) &
            ,safe_adb(NodulNonstElmConc(ielmp),NodulNonstElmConc(ielmp)+NodulNonstElmConc(ielmc)*CPKI)))
!          if(curday==73)write(*,*)NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),CNKI
          CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmn)/CNKI)))
          CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmp)/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        IF(RootNodulElm_pvr(ielmc,L,NZ).GT.ZEROP(NZ))THEN
          FCNPF=AMIN1(1.0_r8 &
            ,SQRT(RootNodulElm_pvr(ielmn,L,NZ)/(RootNodulElm_pvr(ielmc,L,NZ)*NodulerNC_pft(NZ))) &
            ,SQRT(RootNodulElm_pvr(ielmp,L,NZ)/(RootNodulElm_pvr(ielmc,L,NZ)*NodulerPC_pft(NZ))))
        ELSE
          FCNPF=1.0_r8
        ENDIF
        SPNDLI=NodulNonstElmConc(ielmc)/(NodulNonstElmConc(ielmc)+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RespNonst_OUltd=respiration from non-structural C unltd by O2
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDL=bacterial C mass
!     fTgrowRootP_vr=temperature function for root growth
!     FCNPF=N,P constraint to bacterial activity
!     fRootGrowPSISense=growth function of root water potential
!
        RespNonst_OUltd=AZMAX1(AMIN1(RootNodulNonstElm_pvr(ielmc,L,NZ) &
          ,VMXO*RootNodulElm_pvr(ielmc,L,NZ))*FCNPF*fTgrowRootP_vr(L,NZ)*fRootGrowPSISense(1,L))
        CPOOLNX=RootNodulNonstElm_pvr(ielmc,L,NZ)
!
!     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RespNonst_Oltd=respiration from non-structural C ltd by O2
!     RootAutoRO2Limiter_pvr=constraint by O2 consumption on all root processes
!
        RespNonst_Oltd=RespNonst_OUltd*RootAutoRO2Limiter_pvr(1,L,NZ)
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     Rmaint=bacterial maintenance respiration
!     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN6_vr=temperature function for root maintenance respiration
!     WTNDLN=bacterial N mass
!
        Rmaint=AZMAX1(RmSpecPlant*TFN6_vr(L)*RootNodulElm_pvr(ielmn,L,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDLM,RXNDL=difference between non-structural C respn and mntc respn unltd,ltd by O2
!     RGNDLM,RGNDL=growth respiration unlimited by N,P and unltd,ltd by O2
!     RSNDLM,RSNDL=excess maintenance respiration unltd,ltd by O2
!
        RXNDLM=RespNonst_OUltd-Rmaint
        RXNDL=RespNonst_Oltd-Rmaint
        RGNDLM=AZMAX1(RXNDLM)
        RGNDL=AZMAX1(RXNDL)
        RSNDLM=AZMAX1(-RXNDLM)
        RSNDL=AZMAX1(-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDL,WTNDLN=bacterial C,N mass
!     NodulerNC_pft=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     Rauto4Nfix_rvr=respiration for N2 fixation
!     RootN2Fix_pvr,RootN2Fix_pft=layer,total root N2 fixation
!
        RGN2P=AZMAX1(RootNodulElm_pvr(ielmc,L,NZ)*NodulerNC_pft(NZ)-RootNodulElm_pvr(ielmn,L,NZ))/EN2F
        IF(RGNDL.GT.ZEROP(NZ))THEN
          Rauto4Nfix_rvr=RGNDL*RGN2P/(RGNDL+RGN2P)
        ELSE
          Rauto4Nfix_rvr=0._r8
        ENDIF
        RootN2Fix_pvr(L,NZ)=Rauto4Nfix_rvr*EN2F
        RootN2Fix_pft(NZ)=RootN2Fix_pft(NZ)+RootN2Fix_pvr(L,NZ)
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     WTRTD=root C mass
!     CCNDLR=bacteria:root ratio
!     RDNDLX=effect of CCNDLR on bacteria decomposition rate
!     CCNKR=Km for bacterial vs root mass in decomposition
!     SPNDX=specific bacterial decomposition rate at current CCNDLR
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     NoduleElmDecayLoss(ielmc),NoduleElmDecayLoss(ielmn),NoduleElmDecayLoss(ielmp)=bacterial C,N,P loss from decomposition
!     NoduleElmntDecay2Litr(ielmc),NoduleElmntDecay2Litr(ielmn),NoduleElmntDecay2Litr(ielmp)=bacterial C,N,P decomposition to LitrFall
!     NodulElmDecayRecyc(ielmc),NodulElmDecayRecyc(ielmn),NodulElmDecayRecyc(ielmp)=bacterial C,N,P decomposition to recycling
!
        RCCC=RCCZN+CCC*RCCYN
        RCCN=CNC*RCCXN
        RCCP=CPC*RCCQN
        SPNDX=SPNDL*SQRT(fTgrowRootP_vr(L,NZ)*fRootGrowPSISense(1,L))
        DO NE=1,NumPlantChemElms
          NoduleElmDecayLoss(NE)=SPNDX*RootNodulElm_pvr(NE,L,NZ)
        ENDDO

        NoduleElmntDecay2Litr(ielmc)=NoduleElmDecayLoss(ielmc)*(1.0_r8-RCCC)
        NoduleElmntDecay2Litr(ielmn)=NoduleElmDecayLoss(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
        NoduleElmntDecay2Litr(ielmp)=NoduleElmDecayLoss(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
        NodulElmDecayRecyc(ielmc)=NoduleElmDecayLoss(ielmc)-NoduleElmntDecay2Litr(ielmc)
        NodulElmDecayRecyc(ielmn)=NoduleElmDecayLoss(ielmn)-NoduleElmntDecay2Litr(ielmn)
        NodulElmDecayRecyc(ielmp)=NoduleElmDecayLoss(ielmp)-NoduleElmntDecay2Litr(ielmp)
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     NodulUseNonstElm_rvr(ielmc)=total non-structural C used in bacterial growth and growth respiration
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     Rmaint=bacterial maintenance respiration
!     RespNonst_Oltd=respiration from non-structural C
!     NodulElmDecayRecyc(ielmc)=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     Rauto4Nfix_brch=respiration for N2 fixation
!     NoduleBiomCGrowth=bacterial growth
!     NoduGrowthYield_pft=bacterial growth yield
!     NoduleCResp=bacterial respiration for growth and N2 fixation
!     NodulUseNonstElm_rvr(ielmn),NodulUseNonstElm_rvr(ielmp)=nonstructural N,P used in growth
!     CNND,NodulerPC_pft=bacterial N:C,P:C ratio from PFT file
!     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
        NodulUseNonstElm_rvr(ielmc)=AMIN1(RootNodulNonstElm_pvr(ielmc,L,NZ)-AMIN1(Rmaint,RespNonst_Oltd) &
          -Rauto4Nfix_brch+NodulElmDecayRecyc(ielmc),(RGNDL-Rauto4Nfix_brch)/(1.0_r8-NoduGrowthYield_pft(NZ)))
        NoduleBiomCGrowth=NodulUseNonstElm_rvr(ielmc)*NoduGrowthYield_pft(NZ)
        NoduleCResp=Rauto4Nfix_brch+NodulUseNonstElm_rvr(ielmc)*(1.0_r8-NoduGrowthYield_pft(NZ))
        NodulUseNonstElm_rvr(ielmn)=AZMAX1(AMIN1(RootNodulNonstElm_pvr(ielmn,L,NZ),NoduleBiomCGrowth*NodulerNC_pft(NZ)))*NodulNonstElmConc(ielmn)/(NodulNonstElmConc(ielmn)+CZKM)
        NodulUseNonstElm_rvr(ielmp)=AZMAX1(AMIN1(RootNodulNonstElm_pvr(ielmp,L,NZ),NoduleBiomCGrowth*NodulerPC_pft(NZ)))*NodulNonstElmConc(ielmp)/(NodulNonstElmConc(ielmp)+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmp)=bacterial C,N,P loss from senescence
!     NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmp)=bacterial C,N,P senescence to LitrFall
!     NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmp)=bacterial C,N,P senescence to recycling
!
        IF(RSNDL.GT.0.0_r8.AND.RootNodulElm_pvr(ielmc,L,NZ).GT.ZEROP(NZ).AND.RCCC.GT.ZERO)THEN
          NodulELmLoss2Senes(ielmc)=RSNDL/RCCC
          NodulELmLoss2Senes(ielmn)=NodulELmLoss2Senes(ielmc)*RootNodulElm_pvr(ielmn,L,NZ)/RootNodulElm_pvr(ielmc,L,NZ)
          NodulELmLoss2Senes(ielmp)=NodulELmLoss2Senes(ielmc)*RootNodulElm_pvr(ielmp,L,NZ)/RootNodulElm_pvr(ielmc,L,NZ)
          NodulELmSenes2Litr(ielmc)=NodulELmLoss2Senes(ielmc)*(1.0_r8-RCCC)
          NodulELmSenes2Litr(ielmn)=NodulELmLoss2Senes(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
          NodulELmSenes2Litr(ielmp)=NodulELmLoss2Senes(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
          DO NE=1,NumPlantChemElms
            NodulELmSenes2Recyc(NE)=NodulELmLoss2Senes(NE)-NodulELmSenes2Litr(NE)
          ENDDO
        ELSE
          NodulELmLoss2Senes(1:NumPlantChemElms)=0._r8
          NodulELmSenes2Litr(1:NumPlantChemElms)=0._r8
          NodulELmSenes2Recyc(1:NumPlantChemElms)=0._r8
        ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
!     GrossResp_pft,CanopyPlusNoduRespC_pft=total,above-ground PFT respiration
!     Rmaint=bacterial maintenance respiration
!     RespNonst_Oltd=respiration from non-structural C
!     NoduleCResp=bacterial respiration for growth and N2 fixation
!     NodulELmSenes2Recyc(ielmc)=bacterial C senescence to recycling
!     RCO2A_pvr=total root respiration
!     RootRespPotent_pvr,RCO2N_pvr,RCO2A_pvr unlimited by O2,nonstructural C
!
        RCO2TM=AMIN1(Rmaint,RespNonst_OUltd)+RGNDLM+NodulELmSenes2Recyc(ielmc)
        RCO2T=AMIN1(Rmaint,RespNonst_Oltd)+NoduleCResp+NodulELmSenes2Recyc(ielmc)
        RootRespPotent_pvr(ipltroot,L,NZ)=RootRespPotent_pvr(ipltroot,L,NZ)+RCO2TM
        RCO2N_pvr(ipltroot,L,NZ)=RCO2N_pvr(ipltroot,L,NZ)+RCO2T
        RCO2A_pvr(ipltroot,L,NZ)=RCO2A_pvr(ipltroot,L,NZ)-RCO2T
!
!     NODULE LitrFall CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     NoduleElmntDecay2Litr(ielmc),NoduleElmntDecay2Litr(ielmn),NoduleElmntDecay2Litr(ielmp)=bacterial C,N,P decomposition to LitrFall
!     NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmp)=bacterial C,N,P senescence to LitrFall
!
        D6370: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitfalChemElm_pvr(NE,M,k_fine_litr,L,NZ)=LitfalChemElm_pvr(NE,M,k_fine_litr,L,NZ)&
              +CFOPE(NE,iroot,M,NZ)*(NoduleElmntDecay2Litr(NE)+NodulELmSenes2Litr(NE))
          ENDDO
        ENDDO D6370
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     Rmaint=bacterial maintenance respiration
!     RespNonst_Oltd=respiration from non-structural C
!     Rauto4Nfix_rvr=respiration for N2 fixation
!     NodulUseNonstElm_rvr(ielmc)=total non-structural C used in bacterial growth and growth respiration
!     NodulElmDecayRecyc(ielmc),NodulElmDecayRecyc(ielmn),NodulElmDecayRecyc(ielmp)=bacterial C,N,P decomposition to recycling
!     NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmp)=bacterial C,N,P senescence to recycling
!     NodulUseNonstElm_rvr(ielmn),NodulUseNonstElm_rvr(ielmp)=nonstructural N,P used in growth
!     RootN2Fix_pvr=root N2 fixation
!
        RootNodulNonstElm_pvr(ielmc,L,NZ)=RootNodulNonstElm_pvr(ielmc,L,NZ)-AMIN1(Rmaint,RespNonst_Oltd)&
          -Rauto4Nfix_rvr-NodulUseNonstElm_rvr(ielmc)+NodulElmDecayRecyc(ielmc)
        RootNodulNonstElm_pvr(ielmn,L,NZ)=RootNodulNonstElm_pvr(ielmn,L,NZ)-NodulUseNonstElm_rvr(ielmn) &
          +NodulElmDecayRecyc(ielmn)+NodulELmSenes2Recyc(ielmn)+RootN2Fix_pvr(L,NZ)
        RootNodulNonstElm_pvr(ielmp,L,NZ)=RootNodulNonstElm_pvr(ielmp,L,NZ)-NodulUseNonstElm_rvr(ielmp) &
          +NodulElmDecayRecyc(ielmp)+NodulELmSenes2Recyc(ielmp)
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     NoduleBiomCGrowth=bacterial growth
!     NoduleElmDecayLoss(ielmc),NoduleElmDecayLoss(ielmn),NoduleElmDecayLoss(ielmp)=bacterial C,N,P loss from decomposition
!     NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmp)=bacterial C,N,P loss from senescence
!     NodulUseNonstElm_rvr(ielmn),NodulUseNonstElm_rvr(ielmp)=nonstructural N,P used in growth
!
        RootNodulElm_pvr(ielmc,L,NZ)=RootNodulElm_pvr(ielmc,L,NZ)+NoduleBiomCGrowth &
          -NoduleElmDecayLoss(ielmc)-NodulELmLoss2Senes(ielmc)
        RootNodulElm_pvr(ielmn,L,NZ)=RootNodulElm_pvr(ielmn,L,NZ)+NodulUseNonstElm_rvr(ielmn) &
          -NoduleElmDecayLoss(ielmn)-NodulELmLoss2Senes(ielmn)
        RootNodulElm_pvr(ielmp,L,NZ)=RootNodulElm_pvr(ielmp,L,NZ)+NodulUseNonstElm_rvr(ielmp) &
          -NoduleElmDecayLoss(ielmp)-NodulELmLoss2Senes(ielmp)
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOLR,ZPOOLR,PPOOLR=root non-structural C,N,P mass
!     WTRTD=root C mass
!     WTNDL=bacterial C mass
!     NodulBiomCatInfection=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGR=parameter to calculate nonstructural C,N,P exchange
!     CCNDLR=bacteria:root ratio
!     XFRE(ielmc),XFRE(ielmn),XFRE(ielmp)=nonstructural C,N,P transfer
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(RootMycoNonstElm_pvr(ielmc,ipltroot,L,NZ).GT.ZEROP(NZ) &
          .AND. PopuPlantRootC_vr(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
          CCNDLR=RootNodulElm_pvr(ielmc,L,NZ)/PopuPlantRootC_vr(ipltroot,L,NZ)
          WTRTD1=PopuPlantRootC_vr(ipltroot,L,NZ)
          WTNDL1=AMIN1(PopuPlantRootC_vr(ipltroot,L,NZ),AMAX1(NodulBiomCatInfection*AREA3(NU),RootNodulElm_pvr(ielmc,L,NZ)))
          WTRTDT=WTRTD1+WTNDL1
          IF(WTRTDT.GT.ZEROP(NZ))THEN
            FXRNX=FXRN(iPlantNfixType(NZ))/(1.0_r8+CCNDLR/CCNGR)
!    2/(1.0+CCNDLR/(CCNGR*FXRN(iPlantNfixType(NZ))))
            CPOOLD=(RootMycoNonstElm_pvr(ielmc,ipltroot,L,NZ)*WTNDL1-RootNodulNonstElm_pvr(ielmc,L,NZ)*WTRTD1)/WTRTDT
            XFRE(ielmc)=FXRNX*CPOOLD
            RootMycoNonstElm_pvr(ielmc,ipltroot,L,NZ)= RootMycoNonstElm_pvr(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
            RootNodulNonstElm_pvr(ielmc,L,NZ)=RootNodulNonstElm_pvr(ielmc,L,NZ)+XFRE(ielmc)
            CPOOLT= RootMycoNonstElm_pvr(ielmc,ipltroot,L,NZ)+RootNodulNonstElm_pvr(ielmc,L,NZ)
            IF(CPOOLT.GT.ZEROP(NZ))THEN
              ZPOOLD=(RootMycoNonstElm_pvr(ielmn,ipltroot,L,NZ)*RootNodulNonstElm_pvr(ielmc,L,NZ) &
                -RootNodulNonstElm_pvr(ielmn,L,NZ)* RootMycoNonstElm_pvr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(ielmn)=FXRNX*ZPOOLD
              PPOOLD=(RootMycoNonstElm_pvr(ielmp,ipltroot,L,NZ)*RootNodulNonstElm_pvr(ielmc,L,NZ) &
                -RootNodulNonstElm_pvr(ielmp,L,NZ)* RootMycoNonstElm_pvr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(ielmp)=FXRNX*PPOOLD
              DO NE=2,NumPlantChemElms
                RootMycoNonstElm_pvr(NE,ipltroot,L,NZ)=RootMycoNonstElm_pvr(NE,ipltroot,L,NZ)-XFRE(NE)
                RootNodulNonstElm_pvr(NE,L,NZ)=RootNodulNonstElm_pvr(NE,L,NZ)+XFRE(NE)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D5400
  ENDIF
  end associate
  end subroutine RootNodulBiochemistry
end module NoduleBGCMod
