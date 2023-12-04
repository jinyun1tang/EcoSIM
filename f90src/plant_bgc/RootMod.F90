module RootMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod  , only : safe_adb,AZMAX1,AZMIN1
  use EcosimConst
  use GrosubPars
  use PlantMathFuncMod
  use PlantAPIData
  use NoduleBGCMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: RootBGCModel
  contains

  subroutine RootBGCModel(I,J,NZ,BegRemoblize,ICHK1,IDTHRN,NRX,PTRT,TFN6,CNRTW,CPRTW,RootAreaPopu)

  implicit none
  integer, intent(in) :: I,J,NZ,BegRemoblize
  integer, intent(inout) :: ICHK1(jroots,JZ1)
  integer, intent(inout)  :: NRX(jroots,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),CNRTW,CPRTW,RootAreaPopu
  real(r8), intent(in):: PTRT
  integer, intent(out) :: IDTHRN
  real(r8) :: RTVL
  real(r8) :: fRootGrowPsiSense(jroots,JZ1)
  real(r8) :: RootSinkC_vr(jroots,JZ1)
  real(r8) :: Root1stSink_pvr(jroots,JZ1,10)
  real(r8) :: Root2ndSink_pvr(jroots,JZ1,10)
  real(r8) :: RootSinkC(jroots)

  associate(                                &
    NonstructalChemElmnts_pft    =>   plt_biom%NonstructalChemElmnts_pft      , &
    ZEROL    =>   plt_biom%ZEROL      , &
    DLYR3    =>   plt_site%DLYR3      , &
    PlantPopulation_pft       =>   plt_site%PlantPopulation_pft         , &
    ZERO     =>   plt_site%ZERO       , &
    iPlantPhenologyPattern_pft   =>   plt_pheno%iPlantPhenologyPattern_pft    , &
    iPlantRootState_pft    =>   plt_pheno%iPlantRootState_pft     , &
    iPlantShootState_pft    =>   plt_pheno%iPlantShootState_pft     , &
    SeedLengthMean_pft     =>   plt_morph%SeedLengthMean_pft      , &
    RootVH2O_vr   =>   plt_morph%RootVH2O_vr    , &
    RootAreaPerPlant_vr    =>   plt_morph%RootAreaPerPlant_vr     , &
    RootVolume_vr    =>   plt_morph%RootVolume_vr     , &
    RootLenthDensPerPopu_pvr    =>   plt_morph%RootLenthDensPerPopu_pvr     , &
    RootLenPerPopu_pvr    =>   plt_morph%RootLenPerPopu_pvr     , &
    NGTopRootLayer_pft       =>   plt_morph%NGTopRootLayer_pft        , &
    RootPorosity     =>   plt_morph%RootPorosity      , &
    SeedVolumeMean_pft     =>   plt_morph%SeedVolumeMean_pft      , &
    SeedAreaMean_pft     =>   plt_morph%SeedAreaMean_pft      , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft         &
  )
!     ROOT GROWTH
!
  call RootBiochemistry(I,J,NZ,ICHK1,IDTHRN,NRX,TFN6,CNRTW,CPRTW,&
      RootAreaPopu,fRootGrowPsiSense,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)
!
!     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
!     GERMINATION)
!
  RootLenPerPopu_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootLenPerPopu_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedLengthMean_pft(NZ)
  IF(DLYR3(NGTopRootLayer_pft(NZ)).GT.ZERO)THEN
    RootLenthDensPerPopu_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootLenPerPopu_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)/DLYR3(NGTopRootLayer_pft(NZ))
  ELSE
    RootLenthDensPerPopu_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=0._r8
  ENDIF
  RTVL=RootVolume_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RootVH2O_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedVolumeMean_pft(NZ)*PlantPopulation_pft(NZ)
  RootVolume_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootPorosity(ipltroot,NZ)*RTVL
  RootVH2O_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=(1.0_r8-RootPorosity(ipltroot,NZ))*RTVL
  RootAreaPerPlant_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootAreaPerPlant_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+&
    SeedAreaMean_pft(NZ)

  IF(IDTHRN.EQ.NumRootAxes_pft(NZ).OR. &
    (NonstructalChemElmnts_pft(ielmc,NZ).LE.ZEROL(NZ).AND. &
    iPlantPhenologyPattern_pft(NZ).NE.iplt_annual))THEN
    iPlantRootState_pft(NZ)=iDead
    iPlantShootState_pft(NZ)=iDead
  ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
  call RootNoduleBiomchemistry(I,J,NZ,TFN6,fRootGrowPsiSense)

  call NonstructlBiomTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,BegRemoblize)
  end associate
  end subroutine RootBGCModel

!------------------------------------------------------------------------------------------

  subroutine RootBiochemistry(I,J,NZ,ICHK1,IDTHRN,NRX,TFN6,CNRTW,CPRTW,RootAreaPopu,&
    fRootGrowPsiSense,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: ICHK1(jroots,JZ1)
  integer, intent(out) :: IDTHRN
  INTEGER, Intent(inout) :: NRX(jroots,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),CNRTW,CPRTW,RootAreaPopu
  real(r8), intent(out) :: fRootGrowPsiSense(jroots,JZ1)
  REAL(R8), INTENT(OUT)  :: RootSinkC_vr(jroots,JZ1)
  real(r8), INTENT(OUT) :: RootSinkC(jroots)
  real(r8), INTENT(OUT) :: Root1stSink_pvr(jroots,JZ1,10),Root2ndSink_pvr(jroots,JZ1,10)
  integer :: LL,LZ,L1,L,K,lx,M,NR,N,NTG
  real(r8) :: WFNR
  real(r8) :: WFNRG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLD
  real(r8) :: CPOOLX
  real(r8) :: DMRTD
  real(r8) :: FWTRT
  real(r8) :: SoilResit4SecndRootPentration
  real(r8) :: TotSecndRootLen,TotPrimRootLen
  real(r8) :: RTLGX
  real(r8) :: RTLGT
  real(r8) :: RTVL
  real(r8) :: RTAR
  real(r8) :: WTRTTX
  real(r8) :: WTRVCX
  real(r8) :: Root2ndC,Root1stC
  real(r8) :: WTRTLX
  real(r8) :: WTRTTT
  real(r8) :: WTRTT
  real(r8) :: XFRC

!     begin_execution
  associate(                            &
    RootMycoNonstructElmnt_vr   =>   plt_biom%RootMycoNonstructElmnt_vr     , &
    RootStructBiomC_vr   =>   plt_biom%RootStructBiomC_vr     , &
    RootChemElmnts_pft    =>   plt_biom%RootChemElmnts_pft      , &
    ZEROP    =>   plt_biom%ZEROP      , &
    NonstructalChemElmnts_pft    =>   plt_biom%NonstructalChemElmnts_pft      , &
    FWODRE   =>   plt_allom%FWODRE    , &
    RootBiomGrowthYield    =>   plt_allom%RootBiomGrowthYield     , &
    iPlantMorphologyType   =>   plt_pheno%iPlantMorphologyType    , &
    PSIRoot_vr   =>   plt_ew%PSIRoot_vr       , &
    PSIRootTurg_vr    =>   plt_ew%PSIRootTurg_vr        , &
    RootGasLoss_disturb    =>   plt_bgcr%RootGasLoss_disturb    , &
    trcg_rootml     =>   plt_rbgc%trcg_rootml       , &
    trcs_rootml  => plt_rbgc%trcs_rootml, &
    SoilResit4RootPentration     =>   plt_soilchem%SoilResit4RootPentration   , &
    VLSoilPoreMicP     =>   plt_soilchem%VLSoilPoreMicP   , &
    NU       =>   plt_site%NU         , &
    ZERO     =>   plt_site%ZERO       , &
    PlantPopulation_pft       =>   plt_site%PlantPopulation_pft         , &
    ZEROS2   =>   plt_site%ZEROS2     , &
    DLYR3    =>   plt_site%DLYR3      , &
    NL       =>   plt_site%NL         , &
    k_fine_litr=> pltpar%k_fine_litr  , &
    RootVolume_vr    =>   plt_morph%RootVolume_vr     , &
    RootLenthDensPerPopu_pvr    =>   plt_morph%RootLenthDensPerPopu_pvr     , &
    RootAreaPerPlant_vr    =>   plt_morph%RootAreaPerPlant_vr     , &
    RootVH2O_vr   =>   plt_morph%RootVH2O_vr    , &
    Max1stRootRadius1   =>   plt_morph%Max1stRootRadius1    , &
    Max2ndRootRadius1   =>   plt_morph%Max2ndRootRadius1    , &
    Max1stRootRadius   =>   plt_morph%Max1stRootRadius    , &
    PrimRootRadius   =>   plt_morph%PrimRootRadius    , &
    RootLenPerPopu_pvr    =>   plt_morph%RootLenPerPopu_pvr     , &
    AveSecndRootLen    =>   plt_morph%AveSecndRootLen     , &
    SecndRootRadius    =>   plt_morph%SecndRootRadius     , &
    SecndRootXSecArea   =>   plt_morph%SecndRootXSecArea    , &
    Max2ndRootRadius   =>   plt_morph%Max2ndRootRadius    , &
    SecndRootXNum_pvr     =>   plt_morph%SecndRootXNum_pvr      , &
    PrimRootXSecArea   =>   plt_morph%PrimRootXSecArea    , &
    RootPorosity    =>   plt_morph%RootPorosity     , &
    NI       =>   plt_morph%NI        , &
    RootVolPerMassC_pft     =>   plt_morph%RootVolPerMassC_pft      , &
    SeedLengthMean_pft     =>   plt_morph%SeedLengthMean_pft      , &
    MY       =>   plt_morph%MY        , &
    NGTopRootLayer_pft       =>   plt_morph%NGTopRootLayer_pft        , &
    NIXBotRootLayer_pft     =>   plt_morph%NIXBotRootLayer_pft        &
  )

  NIXBotRootLayer_pft(NZ)=NGTopRootLayer_pft(NZ)
  IDTHRN=0
!
  call SummarizeRootSink(NZ,RootAreaPopu,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)
!
!     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
!
  D5010: DO N=1,MY(NZ)
    D5000: DO L=NU,NI(NZ)
!
!     IDENTIFY NEXT LOWER ROOT LAYER
!
!     VLSoilPoreMicP=soil layer volume excluding macropore, rocks
!
      IF(VLSoilPoreMicP(L).GT.ZEROS2)THEN
        D5003: DO LZ=L+1,NL
          IF(VLSoilPoreMicP(LZ).GT.ZEROS2.OR.LZ.EQ.NL)THEN
            L1=LZ
            EXIT
          ENDIF
        ENDDO D5003
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     SoilResit4RootPentration,SoilResit4SecndRootPentration=soil resistance to secondary root penetration (MPa)
!     SecndRootRadius=secondary root radius
!     WFNR=water function for root extension
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     fRootGrowPsiSense,WFNRG=growth,respiration function of root water potential
!     PSIRoot_vr,PSIRootTurg_vr=root total,turgor water potential
!     DMRT=root growth yield
!
        SoilResit4SecndRootPentration=SoilResit4RootPentration(L)*SecndRootRadius(N,L,NZ)/1.0E-03_r8
        WFNR=AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-PSIMin4OrganExtension-SoilResit4SecndRootPentration))
        IF(is_plant_bryophyte(iPlantMorphologyType(NZ)))THEN
          fRootGrowPsiSense(N,L)=EXP(0.05_r8*PSIRoot_vr(N,L,NZ))
          WFNRG=WFNR**0.10_r8
        ELSE
          fRootGrowPsiSense(N,L)=EXP(0.10_r8*PSIRoot_vr(N,L,NZ))
          WFNRG=WFNR**0.25_r8
        ENDIF
        DMRTD=1.0_r8-RootBiomGrowthYield(NZ)
!
!     FOR EACH ROOT AXIS
!
        call GrowRootAxes(N,L,L1,NZ,NRX,fRootGrowPsiSense,ICHK1,WFNR,WFNRG,TFN6,RootAreaPopu,DMRTD,&
          RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,CNRTW,CPRTW,TotPrimRootLen,Root2ndC,Root1stC,TotSecndRootLen)

!
!     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
!     SEASONAL STORAGE POOL IS DEPLETED
!
!     RootStructBiomC_vr,WTRT=total root C mass
!     WTRVC=storage C
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!     CPOOLR=non-structural C mass in root
!
        IF(L.LE.NIXBotRootLayer_pft(NZ))THEN
          IF(RootStructBiomC_vr(N,L,NZ).GT.ZEROP(NZ).AND. &
            RootChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ) &
            .AND.NonstructalChemElmnts_pft(ielmc,NZ).LT.XFRX*RootChemElmnts_pft(ielmc,NZ))THEN
            FWTRT=RootStructBiomC_vr(N,L,NZ)/RootChemElmnts_pft(ielmc,NZ)
            WTRTLX=RootStructBiomC_vr(N,L,NZ)
            WTRTTX=RootChemElmnts_pft(ielmc,NZ)*FWTRT
            WTRTTT=WTRTLX+WTRTTX
            CPOOLX=AZMAX1(RootMycoNonstructElmnt_vr(ielmc,N,L,NZ))
            WTRVCX=AZMAX1(NonstructalChemElmnts_pft(ielmc,NZ)*FWTRT)
            CPOOLD=(WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
            XFRC=AZMIN1(XFRY*CPOOLD)
            RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)+XFRC
            NonstructalChemElmnts_pft(ielmc,NZ)=NonstructalChemElmnts_pft(ielmc,NZ)-XFRC
          ENDIF
        ENDIF
!
!     ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
!     TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
!
!     TotPrimRootLen=total primary root length
!     Root1stC=total primary root C mass
!     TotSecndRootLen=total secondary root length
!     Root2ndC=total secondary root C mass
!     RTLGT=total root length
!     WTRTT=total root C mass
!     FWOOD=C woody fraction in root:0=woody,1=non-woody
!     PP=PFT population
!     RootLenthDensPerPopu_pvr,RootLenPerPopu_pvr=root length density,root length per plant
!     RTVL,RootVH2O_vr,RootVolume_vr=root or myco total,aqueous,gaseous volume
!     RRAD1,SecndRootRadius=primary,secondary root radius
!     RootAreaPerPlant_vr=root surface area per plant
!     AveSecndRootLen=average secondary root length
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!
        IF(N.EQ.ipltroot)THEN
          TotPrimRootLen=TotPrimRootLen*FWODRE(ielmc,k_fine_litr)
          TotSecndRootLen=TotSecndRootLen*FWODRE(ielmc,k_fine_litr)
        ENDIF
        RTLGX=TotPrimRootLen*PlantPopulation_pft(NZ)
        RTLGT=TotSecndRootLen+RTLGX
        WTRTT=Root2ndC+Root1stC
        IF(RTLGT.GT.ZEROP(NZ).AND.WTRTT.GT.ZEROP(NZ).AND.PlantPopulation_pft(NZ).GT.ZEROP(NZ))THEN
          RootLenPerPopu_pvr(N,L,NZ)=RTLGT/PlantPopulation_pft(NZ)
          IF(DLYR3(L).GT.ZERO)THEN
            RootLenthDensPerPopu_pvr(N,L,NZ)=RootLenPerPopu_pvr(N,L,NZ)/DLYR3(L)
          ELSE
            RootLenthDensPerPopu_pvr(N,L,NZ)=0._r8
          ENDIF
          RTVL=AMAX1(PrimRootXSecArea(N,NZ)*RTLGX+SecndRootXSecArea(N,NZ)*TotSecndRootLen &
            ,WTRTT*RootVolPerMassC_pft(N,NZ)*PSIRootTurg_vr(N,L,NZ))
          RootVolume_vr(N,L,NZ)=RootPorosity(N,NZ)*RTVL
          RootVH2O_vr(N,L,NZ)=(1.0_r8-RootPorosity(N,NZ))*RTVL
          !primary roots
          PrimRootRadius(N,L,NZ)=AMAX1(Max1stRootRadius1(N,NZ),&
            (1.0_r8+PSIRoot_vr(N,L,NZ)/EMODR)*Max1stRootRadius(N,NZ))
          !secondary roots
          SecndRootRadius(N,L,NZ)=AMAX1(Max2ndRootRadius1(N,NZ),&
            (1.0_r8+PSIRoot_vr(N,L,NZ)/EMODR)*Max2ndRootRadius(N,NZ))
          RTAR=TwoPiCON*(PrimRootRadius(N,L,NZ)*RTLGX+SecndRootRadius(N,L,NZ)*TotSecndRootLen)
          IF(SecndRootXNum_pvr(N,L,NZ).GT.ZEROP(NZ))THEN
            AveSecndRootLen(N,L,NZ)=AMAX1(MinAve2ndRootLen,TotSecndRootLen/SecndRootXNum_pvr(N,L,NZ))
          ELSE
            AveSecndRootLen(N,L,NZ)=MinAve2ndRootLen
          ENDIF
          RootAreaPerPlant_vr(N,L,NZ)=RTAR/PlantPopulation_pft(NZ)
!     IF(N.EQ.ipltroot)THEN
!     RootAreaPerPlant_vr(N,L,NZ)=RootAreaPerPlant_vr(N,L,NZ)*MinAve2ndRootLen/AveSecndRootLen(N,L,NZ)
!     ENDIF
        ELSE
          RootLenPerPopu_pvr(N,L,NZ)=0._r8
          RootLenthDensPerPopu_pvr(N,L,NZ)=0._r8
          RootVolume_vr(N,L,NZ)=0._r8
          RootVH2O_vr(N,L,NZ)=0._r8
          PrimRootRadius(N,L,NZ)=Max1stRootRadius(N,NZ)
          SecndRootRadius(N,L,NZ)=Max2ndRootRadius(N,NZ)
          RootAreaPerPlant_vr(N,L,NZ)=0._r8
          AveSecndRootLen(N,L,NZ)=MinAve2ndRootLen
          DO NTG=idg_beg,idg_end-1
            RootGasLoss_disturb(NTG,NZ)=RootGasLoss_disturb(NTG,NZ)-(trcg_rootml(NTG,N,L,NZ)+trcs_rootml(NTG,N,L,NZ))
          ENDDO
          trcg_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
          trcs_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
        ENDIF
      ENDIF
    ENDDO D5000
  ENDDO D5010
  end associate
  end subroutine RootBiochemistry

!------------------------------------------------------------------------------------------

  subroutine GrowRootAxes(N,L,L1,NZ,NRX,fRootGrowPsiSense,ICHK1,WFNR,WFNRG,TFN6,RootAreaPopu,DMRTD,RootSinkC_vr,&
    Root1stSink_pvr,Root2ndSink_pvr,CNRTW,CPRTW,TotSecndRootLen,TotPrimRootLen,Root2ndC,Root1stC)
  implicit none
  INTEGER, INTENT(IN) :: N,L,L1,NZ
  integer, intent(inout) :: NRX(2,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),RootAreaPopu
  real(r8), intent(in) :: fRootGrowPsiSense(2,JZ1)
  real(r8), intent(in) :: DMRTD
  REAL(R8), INTENT(IN) :: RootSinkC_vr(2,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(2,JZ1,10),Root2ndSink_pvr(2,JZ1,10),CNRTW,CPRTW
  integer, intent(inout) :: ICHK1(2,JZ1)
  real(r8), intent(inout):: WFNR,WFNRG
  real(r8), intent(out) :: TotSecndRootLen,TotPrimRootLen,Root2ndC,Root1stC
  real(r8) :: CNRDA,CNRDM
  real(r8) :: CNPG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: FSNC1,FSNC2
  real(r8) :: ZADD1,ZADD1M,ZADD2,ZADD2M
  real(r8) :: ZPOOLB
  real(r8) :: CGRORM
  real(r8) :: CGROR
  real(r8) :: DMRTR
  real(r8) :: FNP
  real(r8) :: FRTN
  real(r8) :: FRCO2
  real(r8) :: FSNCM
  real(r8) :: FSNCP
  real(r8) :: GRTWGM
  real(r8) :: Root1stExtension
  real(r8) :: RootCYieldO2ltd
  real(r8) :: RootNLigthSatCarboxyRate_nodewthElmnt(NumOfPlantChemElmnts)
  real(r8) :: GRTWTM
  real(r8) :: PPOOLB
  real(r8) :: PADD2,PADD1
  real(r8) :: RCO2X,RCO2Y
  real(r8) :: RCO2G,RCO2T
  real(r8) :: RCO2XM
  real(r8) :: RCO2YM
  real(r8) :: RCO2GM
  real(r8) :: RCO2TM
  real(r8) :: RMNCR,RCO2RM,RCO2R
  real(r8) :: RootChemElmntRemob(NumOfPlantChemElmnts)
  real(r8) :: RTN2X,RTN2Y
  real(r8) :: RTDP1X,SoilResit4PrimRootPentration
  REAL(R8) :: SNCR,SNCRM
  real(r8) :: TFRCO2
  real(r8) :: RCCC,RCCN,RCCP
  integer :: NR,M,LL,LX,NE

!begin_execution
  associate(                          &
    Root1stChemElmnt  =>  plt_biom%Root1stChemElmnt     , &
    Root2ndStructChemElmnt_pvr  =>  plt_biom%Root2ndStructChemElmnt_pvr     , &
    RootNonstructElementConcpft_vr  =>  plt_biom%RootNonstructElementConcpft_vr     , &
    Root1stStructChemElmnt_pvr  =>  plt_biom%Root1stStructChemElmnt_pvr     , &
    RootMycoNonstructElmnt_vr  =>  plt_biom%RootMycoNonstructElmnt_vr     , &
    RootProteinC_pvr   =>  plt_biom%RootProteinC_pvr      , &
    RootStructBiomC_vr  =>  plt_biom%RootStructBiomC_vr     , &
    ZEROP   =>  plt_biom%ZEROP      , &
    CumSoilThickness  =>  plt_site%CumSoilThickness     , &
    RCO2A   =>  plt_rbgc%RCO2A      , &
    RCO2N   =>  plt_rbgc%RCO2N      , &
    RootRespPotential_vr   =>  plt_rbgc%RootRespPotential_vr      , &
    RootAutoRO2Limiter_pvr    =>  plt_rbgc%RootAutoRO2Limiter_pvr       , &
    LitterFallChemElmnt_pftvr    =>  plt_bgcr%LitterFallChemElmnt_pftvr       , &
    rCNNonstructRemob_pft   =>  plt_allom%rCNNonstructRemob_pft     , &
    rCPNonstructRemob_pft    =>  plt_allom%rCPNonstructRemob_pft      , &
    FWODRE  =>  plt_allom%FWODRE    , &
    CNRTS   =>  plt_allom%CNRTS     , &
    CPRTS   =>  plt_allom%CPRTS     , &
    RootBiomGrowthYield    =>  plt_allom%RootBiomGrowthYield      , &
    k_woody_litr=> pltpar%k_woody_litr,&
    k_fine_litr=> pltpar%k_fine_litr, &
    icwood  =>  pltpar%icwood       , &
    iroot   =>  pltpar%iroot        , &
    inonstruct=>  pltpar%inonstruct     , &
    iPlantMorphologyType  =>  plt_pheno%iPlantMorphologyType    , &
    iPlantPhenologyType_pft  =>  plt_pheno%iPlantPhenologyType_pft    , &
    fTgrowRootP    =>  plt_pheno%fTgrowRootP      , &
    iPlantCalendar_brch  =>  plt_pheno%iPlantCalendar_brch    , &
    SoiBulkDensity    =>  plt_soilchem%SoiBulkDensity   , &
    CFOPE   =>  plt_soilchem%CFOPE  , &
    SoilResit4RootPentration    =>  plt_soilchem%SoilResit4RootPentration   , &
    DLYR3   =>  plt_site%DLYR3      , &
    ZERO    =>  plt_site%ZERO       , &
    MaxNumRootLays      =>  plt_site%MaxNumRootLays         , &
    PSIRootTurg_vr   =>  plt_ew%PSIRootTurg_vr        , &
    SecndRootXNum_pvr    =>  plt_morph%SecndRootXNum_pvr      , &
    MaxSeedCMass    =>  plt_morph%MaxSeedCMass      , &
    PrimRootDepth  =>  plt_morph%PrimRootDepth    , &
    PrimRootXNumL   =>  plt_morph%PrimRootXNumL     , &
    NGTopRootLayer_pft     =>  plt_morph%NGTopRootLayer_pft       , &
    NIXBotRootLayer_pft    =>  plt_morph%NIXBotRootLayer_pft      , &
    SecndRootXNum_rpvr   =>  plt_morph%SecndRootXNum_rpvr     , &
    NumRootAxes_pft    =>  plt_morph%NumRootAxes_pft      , &
    PrimRootRadius  =>  plt_morph%PrimRootRadius    , &
    RootBranchFreq_pft    =>  plt_morph%RootBranchFreq_pft      , &
    SeedinDepth   =>  plt_morph%SeedinDepth     , &
    PrimRootLen   =>  plt_morph%PrimRootLen     , &
    SecndRootSpecLen  =>  plt_morph%SecndRootSpecLen    , &
    SecndRootLen   =>  plt_morph%SecndRootLen     , &
    NIXBotRootLayer_rpft    =>  plt_morph%NIXBotRootLayer_rpft      , &
    NumOfMainBranch_pft     =>  plt_morph%NumOfMainBranch_pft       , &
    C4PhotosynDowreg_brch   =>  plt_photo%C4PhotosynDowreg_brch       &
  )
  TotSecndRootLen=0._r8
  TotPrimRootLen=0._r8
  Root2ndC=0._r8
  Root1stC=0._r8
  D5050: DO NR=1,NumRootAxes_pft(NZ)
!
!     SECONDARY ROOT EXTENSION
!
    IF(L.LE.NIXBotRootLayer_rpft(NR,NZ).AND.NRX(N,NR).EQ.0)THEN
!
!     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
!     TO CURRENT AXIS
!
!     Root2ndSink_pvr=total secondary root sink strength
!     RootSinkC_vr=total root sink strength
!     FRTN=fraction of secondary root sink strength in axis
!
      IF(RootSinkC_vr(N,L).GT.ZEROP(NZ))THEN
        FRTN=Root2ndSink_pvr(N,L,NR)/RootSinkC_vr(N,L)
      ELSE
        FRTN=1.0_r8
      ENDIF
!
!     N,P CONSTRAINT ON SECONDARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
      IF(RootNonstructElementConcpft_vr(ielmc,N,L,NZ).GT.ZERO)THEN
        CNPG=AMIN1(RootNonstructElementConcpft_vr(ielmn,N,L,NZ)/(RootNonstructElementConcpft_vr(ielmn,N,L,NZ) &
          +RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CNKI),RootNonstructElementConcpft_vr(ielmp,N,L,NZ) &
          /(RootNonstructElementConcpft_vr(ielmp,N,L,NZ)+RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CPKI))
      ELSE
        CNPG=1.0_r8
      ENDIF
!
!     SECONDARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT2N=secondary root N mass
!     TFN6=temperature function for root maintenance respiration
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iPlantPhenologyType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     fRootGrowPsiSense=growth function of root water potential
!
      RMNCR=AZMAX1(RMPLT*Root2ndStructChemElmnt_pvr(ielmn,N,L,NR,NZ))*TFN6(L)
      IF(is_plant_bryophyte(iPlantMorphologyType(NZ)).OR. &
        iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
        RMNCR=RMNCR*fRootGrowPsiSense(N,L)
      ENDIF
!
!     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of secondary root sink strength in axis
!     CPOOL=non-structural C mass
!     fTgrowRootP=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!     fRootGrowPsiSense=growth function of root water potential
!
      RCO2RM=AZMAX1(VMXC*FRTN* RootMycoNonstructElmnt_vr(ielmc,N,L,NZ) &
        *fTgrowRootP(L,NZ))*CNPG*C4PhotosynDowreg_brch(NumOfMainBranch_pft(NZ),NZ)*fRootGrowPsiSense(N,L)
!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     RootAutoRO2Limiter_pvr=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RCO2R=RCO2RM*RootAutoRO2Limiter_pvr(N,L,NZ)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AZMAX1(RCO2XM)*WFNRG
      RCO2Y=AZMAX1(RCO2X)*WFNRG
!
!     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AZMAX1(RootMycoNonstructElmnt_vr(ielmn,N,L,NZ))
      PPOOLB=AZMAX1(RootMycoNonstructElmnt_vr(ielmp,N,L,NZ))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ),PPOOLB*DMRTR/CPRTS(NZ))
      IF(RCO2YM.GT.0.0_r8)THEN
        RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
        RCO2GM=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
        RCO2G=AMIN1(RCO2Y,FNP*RootAutoRO2Limiter_pvr(N,L,NZ))
      ELSE
        RCO2G=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,RootCYieldO2ltd=root C growth unltd,ltd by O2
!     RootBiomGrowthYield=root growth yield
!     ZADD2M,ZADD2,PADD2=nonstructural N,P unlimited,limited by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
      CGRORM=RCO2GM/DMRTD
      CGROR=RCO2G/DMRTD
      GRTWGM=CGRORM*RootBiomGrowthYield(NZ)
      RootCYieldO2ltd=CGROR*RootBiomGrowthYield(NZ)
      ZADD2M=AZMAX1(GRTWGM*CNRTW)
      ZADD2=AZMAX1(AMIN1(FRTN* RootMycoNonstructElmnt_vr(ielmn,N,L,NZ),RootCYieldO2ltd*CNRTW))
      PADD2=AZMAX1(AMIN1(FRTN* RootMycoNonstructElmnt_vr(ielmp,N,L,NZ),RootCYieldO2ltd*CPRTW))
      CNRDM=AZMAX1(1.70_r8*ZADD2M)
      CNRDA=AZMAX1(1.70_r8*ZADD2)
!
!     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
      IF(iPlantCalendar_brch(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).NE.0.AND.&
        RootNonstructElementConcpft_vr(ielmc,N,L,NZ).GT.ZERO)THEN
        CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElementConcpft_vr(ielmn,N,L,NZ),RootNonstructElementConcpft_vr(ielmn,N,L,NZ) &
          +RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CNKI) &
          ,safe_adb(RootNonstructElementConcpft_vr(ielmp,N,L,NZ),RootNonstructElementConcpft_vr(ielmp,N,L,NZ) &
          +RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CPKI)))
        CNC=AZMAX1(AMIN1(1.0_r8 &
          ,safe_adb(RootNonstructElementConcpft_vr(ielmc,N,L,NZ),RootNonstructElementConcpft_vr(ielmc,N,L,NZ)+ &
          RootNonstructElementConcpft_vr(ielmn,N,L,NZ)/CNKI)))
        CPC=AZMAX1(AMIN1(1.0_r8 &
          ,safe_adb(RootNonstructElementConcpft_vr(ielmc,N,L,NZ),RootNonstructElementConcpft_vr(ielmc,N,L,NZ)+ &
          RootNonstructElementConcpft_vr(ielmp,N,L,NZ)/CPKI)))
      ELSE
        CCC=0._r8
        CNC=0._r8
        CPC=0._r8
      ENDIF
      RCCC=RCCZR+CCC*RCCYR
      RCCN=CNC*RCCXR
      RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P FROM SECONDARY ROOT DURING
!     REMOBILIZATION DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RootAutoRO2Limiter_pvr=constraint by O2 consumption on all root processes
!     RootChemElmntRemob(ielmc),RootChemElmntRemob(ielmn),RootChemElmntRemob(ielmp)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC2=fraction of secondary root C to be remobilized
!
      IF(-RCO2XM.GT.0.0_r8)THEN
        IF(-RCO2XM.LT.Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)*RCCC)THEN
          SNCRM=-RCO2XM
        ELSE
          SNCRM=AZMAX1(Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)*RCCC)
        ENDIF
      ELSE
        SNCRM=0._r8
      ENDIF
      IF(-RCO2X.GT.0.0_r8)THEN
        IF(-RCO2X.LT.Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)*RCCC)THEN
          SNCR=-RCO2X
        ELSE
          SNCR=AZMAX1(Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)*RCCC)*RootAutoRO2Limiter_pvr(N,L,NZ)
        ENDIF
      ELSE
        SNCR=0._r8
      ENDIF
      IF(SNCR.GT.0.0.AND.Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ).GT.ZEROP(NZ))THEN
        RootChemElmntRemob(ielmc)=RCCC*Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
        RootChemElmntRemob(ielmn)=Root2ndStructChemElmnt_pvr(ielmn,N,L,NR,NZ)*(RCCN+(1.0_r8-RCCN)* &
          RootChemElmntRemob(ielmc)/Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ))
        RootChemElmntRemob(ielmp)=Root2ndStructChemElmnt_pvr(ielmp,N,L,NR,NZ)*(RCCP+(1.0_r8-RCCP)* &
          RootChemElmntRemob(ielmc)/Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ))
        IF(RootChemElmntRemob(ielmc).GT.ZEROP(NZ))THEN
          FSNC2=AZMAX1(AMIN1(1.0_r8,SNCR/RootChemElmntRemob(ielmc)))
        ELSE
          FSNC2=1.0_r8
        ENDIF
      ELSE
        RootChemElmntRemob(1:NumOfPlantChemElmnts)=0._r8
        FSNC2=0._r8
      ENDIF
!
!     SECONDARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC2=fraction of secondary root C to be remobilized
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RootChemElmntRemob(ielmc),RootChemElmntRemob(ielmn),RootChemElmntRemob(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      DO NE=1,NumOfPlantChemElmnts
        D6350: DO M=1,jsken
          LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
            *FSNC2*(Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)-RootChemElmntRemob(NE))*FWODRE(NE,k_woody_litr)

          LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
            *FSNC2*(Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)-RootChemElmntRemob(NE))*FWODRE(NE,k_fine_litr)
        ENDDO D6350
      ENDDO
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC2=fraction of secondary root C to be remobilized
!     RootChemElmntRemob(ielmc),RootChemElmntRemob(ielmn),RootChemElmntRemob(ielmp)=remobilization of C,N,P from senescing root
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
       RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)-AMIN1(RMNCR,RCO2R) &
        -CGROR-CNRDA-SNCR+FSNC2*RootChemElmntRemob(ielmc)
       RootMycoNonstructElmnt_vr(ielmn,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmn,N,L,NZ)-ZADD2+FSNC2*RootChemElmntRemob(ielmn)
       RootMycoNonstructElmnt_vr(ielmp,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmp,N,L,NZ)-PADD2+FSNC2*RootChemElmntRemob(ielmp)
!
!     TOTAL SECONDARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RootRespPotential_vr,RCO2N=RCO2A unltd by O2,nonstructural C
!
      RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
      RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
      RootRespPotential_vr(N,L,NZ)=RootRespPotential_vr(N,L,NZ)+RCO2TM
      RCO2N(N,L,NZ)=RCO2N(N,L,NZ)+RCO2T
      RCO2A(N,L,NZ)=RCO2A(N,L,NZ)-RCO2T
!
!     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     Root1stExtension=secondary root length extension
!     RootCYieldO2ltd=secondary root C growth ltd by O2
!     SecndRootSpecLen=specific secondary root length from startq.f
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     FSNC2=fraction of secondary root C to be remobilized
!     SecndRootLen=secondary root length
!     RootNLigthSatCarboxyRate_nodewthElmnt(ielmc),RootNLigthSatCarboxyRate_nodewthElmnt(ielmn),RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=net root C,N,P growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      Root1stExtension=RootCYieldO2ltd*SecndRootSpecLen(N,NZ)*WFNR*FWODRE(ielmc,k_fine_litr) &
        -FSNC2*SecndRootLen(N,L,NR,NZ)
      RootNLigthSatCarboxyRate_nodewthElmnt(ielmc)=RootCYieldO2ltd-FSNC2*Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
      RootNLigthSatCarboxyRate_nodewthElmnt(ielmn)=ZADD2-FSNC2*Root2ndStructChemElmnt_pvr(ielmn,N,L,NR,NZ)
      RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=PADD2-FSNC2*Root2ndStructChemElmnt_pvr(ielmp,N,L,NR,NZ)
!
!     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
!     AND AXIS NUMBER
!
!     SecndRootLen=secondary root length
!     Root1stExtension=secondary root length extension
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RootNLigthSatCarboxyRate_nodewthElmnt(ielmc),RootNLigthSatCarboxyRate_nodewthElmnt(ielmn),RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=net root C,N,P growth
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!     RootBranchFreq_pft=root branching frequency from PFT file
!     SecndRootXNum_rpvr,SecndRootXNum_pvr=number of secondary root axes
!
      SecndRootLen(N,L,NR,NZ)=SecndRootLen(N,L,NR,NZ)+Root1stExtension
      DO NE=1,NumOfPlantChemElmnts
        Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)=Root2ndStructChemElmnt_pvr(NE,N,L,NR,NZ)+RootNLigthSatCarboxyRate_nodewthElmnt(NE)
      ENDDO
      RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)+AMIN1(rCNNonstructRemob_pft(NZ)*Root2ndStructChemElmnt_pvr(ielmn,N,L,NR,NZ) &
        ,rCPNonstructRemob_pft(NZ)*Root2ndStructChemElmnt_pvr(ielmp,N,L,NR,NZ))
      TotSecndRootLen=TotSecndRootLen+SecndRootLen(N,L,NR,NZ)
      Root2ndC=Root2ndC+Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
      RTN2X=RootBranchFreq_pft(NZ)*RootAreaPopu
      RTN2Y=RootBranchFreq_pft(NZ)*RTN2X
      SecndRootXNum_rpvr(N,L,NR,NZ)=(RTN2X+RTN2Y)*DLYR3(L)
      SecndRootXNum_pvr(N,L,NZ)=SecndRootXNum_pvr(N,L,NZ)+SecndRootXNum_rpvr(N,L,NR,NZ)
!
!     PRIMARY ROOT EXTENSION
!
!     SoiBulkDensity=soil bulk density
!     RTDP1,RTDP1X=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     ICHKL=flag for identifying layer with primary root tip
!     RTN1=number of primary root axes
!     RootAreaPopu=multiplier for number of primary root axes
!
      IF(N.EQ.ipltroot)THEN
        IF(SoiBulkDensity(L).GT.ZERO)THEN
          RTDP1X=PrimRootDepth(N,NR,NZ)-CumSoilThickness(0)
        ELSE
          RTDP1X=PrimRootDepth(N,NR,NZ)
        ENDIF
        IF(RTDP1X.GT.CumSoilThickness(L-1).AND.ICHK1(N,NR).EQ.0)THEN
            PrimRootXNumL(N,L,NZ)=PrimRootXNumL(N,L,NZ)+RootAreaPopu
            IF(RTDP1X.LE.CumSoilThickness(L).OR.L.EQ.MaxNumRootLays)THEN
              ICHK1(N,NR)=1
!
!     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
!     ATTRIBUTED TO CURRENT AXIS
!
!     Root1stSink_pvr=primary root sink strength
!     RootSinkC_vr=total root sink strength
!     FRTN=fraction of primary root sink strength in axis
!
              IF(RootSinkC_vr(N,L).GT.ZEROP(NZ))THEN
                FRTN=Root1stSink_pvr(N,L,NR)/RootSinkC_vr(N,L)
              ELSE
                FRTN=1.0_r8
              ENDIF
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     SoilResit4RootPentration,SoilResit4PrimRootPentration=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
!     WFNR=water function for root extension
!     WFNRG=respiration function of root water potential
!
              SoilResit4PrimRootPentration=SoilResit4RootPentration(L)*PrimRootRadius(N,L,NZ)/1.0E-03_r8
              WFNR=AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-PSIMin4OrganExtension-SoilResit4PrimRootPentration))
              IF(is_plant_bryophyte(iPlantMorphologyType(NZ)))THEN
                WFNRG=WFNR**0.10_r8
              ELSE
                WFNRG=WFNR**0.25_r8
              ENDIF
!
!     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
              IF(RootNonstructElementConcpft_vr(ielmc,N,L,NZ).GT.ZERO)THEN
                CNPG=AMIN1(RootNonstructElementConcpft_vr(ielmn,N,L,NZ)/(RootNonstructElementConcpft_vr(ielmn,N,L,NZ) &
                  +RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CNKI),RootNonstructElementConcpft_vr(ielmp,N,L,NZ) &
                  /(RootNonstructElementConcpft_vr(ielmp,N,L,NZ)+RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CPKI))
              ELSE
                CNPG=1.0_r8
              ENDIF
!
!     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT1N=primary root N mass
!     TFN6=temperature function for root maintenance respiration
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iPlantPhenologyType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     fRootGrowPsiSense=growth function of root water potential
!
              RMNCR=AZMAX1(RMPLT*Root1stChemElmnt(ielmn,N,NR,NZ))*TFN6(L)
              IF(is_plant_bryophyte(iPlantMorphologyType(NZ)).OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
                RMNCR=RMNCR*fRootGrowPsiSense(N,L)
              ENDIF
!
!     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of primary root sink strength in axis
!     CPOOL=non-structural C mass
!     fTgrowRootP=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!     fRootGrowPsiSense=growth function of root water potential
!
              RCO2RM=AZMAX1(VMXC*FRTN* RootMycoNonstructElmnt_vr(ielmc,N,L,NZ) &
                *fTgrowRootP(L,NZ))*CNPG*C4PhotosynDowreg_brch(NumOfMainBranch_pft(NZ),NZ)*fRootGrowPsiSense(N,L)
              IF(RTDP1X.GE.CumSoilThickness(MaxNumRootLays))THEN
                RCO2RM=AMIN1(RMNCR,RCO2RM)
              ENDIF
!
!     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     RootAutoRO2Limiter_pvr=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
              RCO2R=RCO2RM*RootAutoRO2Limiter_pvr(N,L,NZ)
              RCO2XM=RCO2RM-RMNCR
              RCO2X=RCO2R-RMNCR
              RCO2YM=AZMAX1(RCO2XM)*WFNRG
              RCO2Y=AZMAX1(RCO2X)*WFNRG
!
!     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
              DMRTR=DMRTD*FRTN
              ZPOOLB=AZMAX1(RootMycoNonstructElmnt_vr(ielmn,N,L,NZ))
              PPOOLB=AZMAX1(RootMycoNonstructElmnt_vr(ielmp,N,L,NZ))
              FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ),PPOOLB*DMRTR/CPRTS(NZ))
              IF(RCO2YM.GT.0.0_r8)THEN
                RCO2GM=AMIN1(RCO2YM,FNP)
              ELSE
                RCO2GM=0._r8
              ENDIF
              IF(RCO2Y.GT.0.0_r8)THEN
                RCO2G=AMIN1(RCO2Y,FNP*RootAutoRO2Limiter_pvr(N,L,NZ))
              ELSE
                RCO2G=0._r8
              ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,RootCYieldO2ltd=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD1M,ZADD1,PADD1=nonstructural N,P unltd,ltd by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
              CGRORM=RCO2GM/DMRTD
              CGROR=RCO2G/DMRTD
              GRTWGM=CGRORM*RootBiomGrowthYield(NZ)
              RootCYieldO2ltd=CGROR*RootBiomGrowthYield(NZ)
              ZADD1M=AZMAX1(GRTWGM*CNRTW)
              ZADD1=AZMAX1(AMIN1(FRTN* RootMycoNonstructElmnt_vr(ielmn,N,L,NZ),RootCYieldO2ltd*CNRTW))
              PADD1=AZMAX1(AMIN1(FRTN* RootMycoNonstructElmnt_vr(ielmp,N,L,NZ),RootCYieldO2ltd*CPRTW))
              CNRDM=AZMAX1(1.70_r8*ZADD1M)
              CNRDA=AZMAX1(1.70_r8*ZADD1)

              call PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC1=fraction of primary root C to be remobilized
!     RootChemElmntRemob(ielmc),RootChemElmntRemob(ielmn),RootChemElmntRemob(ielmp)=remobilization of C,N,P from senescing root
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!
               RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)-AMIN1(RMNCR,RCO2R)-CGROR-CNRDA-SNCR+FSNC1*RootChemElmntRemob(ielmc)
               RootMycoNonstructElmnt_vr(ielmn,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmn,N,L,NZ)-ZADD1+FSNC1*RootChemElmntRemob(ielmn)
               RootMycoNonstructElmnt_vr(ielmp,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmp,N,L,NZ)-PADD1+FSNC1*RootChemElmntRemob(ielmp)
!
!     TOTAL PRIMARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RootRespPotential_vr,RCO2N=RCO2A unltd by O2,nonstructural C
!
              RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
              RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
!
!     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
!     THROUGH WHICH PRIMARY ROOTS GROW
!
!     RTDP1=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     PrimRootLen=primary root length
!     SeedinDepth=seeding depth
!     FRCO2=fraction of primary root respiration attributed to layer
!     RCO2A=total root respiration
!     RootRespPotential_vr,RCO2N=RCO2A unltd by O2,nonstructural C
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!
              IF(PrimRootDepth(N,NR,NZ).GT.CumSoilThickness(NGTopRootLayer_pft(NZ)))THEN
                TFRCO2=0._r8
                D5100: DO LL=NGTopRootLayer_pft(NZ),NIXBotRootLayer_rpft(NR,NZ)
                  IF(LL.LT.NIXBotRootLayer_rpft(NR,NZ))THEN
                    FRCO2=AMIN1(1.0_r8,PrimRootLen(N,LL,NR,NZ)/(PrimRootDepth(N,NR,NZ)-SeedinDepth(NZ)))
                  ELSE
                    FRCO2=1.0_r8-TFRCO2
                  ENDIF
                  TFRCO2=TFRCO2+FRCO2
                  RootRespPotential_vr(N,LL,NZ)=RootRespPotential_vr(N,LL,NZ)+RCO2TM*FRCO2
                  RCO2N(N,LL,NZ)=RCO2N(N,LL,NZ)+RCO2T*FRCO2
                  RCO2A(N,LL,NZ)=RCO2A(N,LL,NZ)-RCO2T*FRCO2
                ENDDO D5100
              ELSE
                RootRespPotential_vr(N,L,NZ)=RootRespPotential_vr(N,L,NZ)+RCO2TM
                RCO2N(N,L,NZ)=RCO2N(N,L,NZ)+RCO2T
                RCO2A(N,L,NZ)=RCO2A(N,L,NZ)-RCO2T
              ENDIF
!
!     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
!     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
!     HAVE DISAPPEARED
!
!     RootCYieldO2ltd=primary root C growth ltd by O2
!     RootNLigthSatCarboxyRate_nodewthElmnt(ielmc),RootNLigthSatCarboxyRate_nodewthElmnt(ielmn),RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=net primary root C,N,P growth
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     SecndRootLen=secondary root length
!
              RootNLigthSatCarboxyRate_nodewthElmnt(ielmc)=RootCYieldO2ltd-FSNC1*Root1stChemElmnt(ielmc,N,NR,NZ)
              RootNLigthSatCarboxyRate_nodewthElmnt(ielmn)=ZADD1-FSNC1*Root1stChemElmnt(ielmn,N,NR,NZ)
              RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=PADD1-FSNC1*Root1stChemElmnt(ielmp,N,NR,NZ)
              IF(RootNLigthSatCarboxyRate_nodewthElmnt(ielmc).LT.0.0_r8)THEN
                LX=MAX(1,L-1)
                D5105: DO LL=L,LX,-1
                  GRTWTM=RootNLigthSatCarboxyRate_nodewthElmnt(ielmc)
                  DO NE=1,NumOfPlantChemElmnts
                    IF(RootNLigthSatCarboxyRate_nodewthElmnt(NE).LT.0.0_r8)THEN
                      IF(RootNLigthSatCarboxyRate_nodewthElmnt(NE).GT.-Root2ndStructChemElmnt_pvr(NE,N,LL,NR,NZ))THEN
                        if(NE==ielmc)SecndRootLen(N,LL,NR,NZ)=SecndRootLen(N,LL,NR,NZ)+RootNLigthSatCarboxyRate_nodewthElmnt(NE) &
                          *SecndRootLen(N,LL,NR,NZ)/Root2ndStructChemElmnt_pvr(NE,N,LL,NR,NZ)
                        Root2ndStructChemElmnt_pvr(NE,N,LL,NR,NZ)=Root2ndStructChemElmnt_pvr(NE,N,LL,NR,NZ)+RootNLigthSatCarboxyRate_nodewthElmnt(NE)
                        RootNLigthSatCarboxyRate_nodewthElmnt(NE)=0._r8
                      ELSE
                        if(NE==ielmc)SecndRootLen(N,LL,NR,NZ)=0._r8
                        RootNLigthSatCarboxyRate_nodewthElmnt(NE)=RootNLigthSatCarboxyRate_nodewthElmnt(NE)+Root2ndStructChemElmnt_pvr(NE,N,LL,NR,NZ)
                        Root2ndStructChemElmnt_pvr(NE,N,LL,NR,NZ)=0._r8
                      ENDIF
                    ENDIF
                  ENDDO
!
!     CONCURRENT LOSS OF MYCORRHIZAE AND NODULES WITH LOSS
!     OF SECONDARY ROOTS
!
!     GRTWTM=negative primary root C growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     FSNCM,FSNCP=fraction of mycorrhizal structural,nonstructural C to be remobilized
!     WTRTL=active root C mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTRT2,WTRT2N,WTRT2P=mycorrhizal C,N,P mass
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     SecndRootLen=mycorrhizal length
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
!
                  IF(GRTWTM.LT.0.0_r8)THEN
                    IF(Root2ndStructChemElmnt_pvr(ielmc,ipltroot,LL,NR,NZ).GT.ZEROP(NZ))THEN
                      FSNCM=AMIN1(1.0_r8,ABS(GRTWTM)/Root2ndStructChemElmnt_pvr(ielmc,ipltroot,LL,NR,NZ))
                    ELSE
                      FSNCM=1.0_r8
                    ENDIF
                    IF(RootStructBiomC_vr(ipltroot,LL,NZ).GT.ZEROP(NZ))THEN
                      FSNCP=AMIN1(1.0_r8,ABS(GRTWTM)/RootStructBiomC_vr(ipltroot,LL,NZ))
                    ELSE
                      FSNCP=1.0_r8
                    ENDIF

                    DO NE=1,NumOfPlantChemElmnts
                      D6450: DO M=1,jsken
                        LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,LL,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,LL,NZ)+CFOPE(NE,icwood,M,NZ) &
                          *FSNCM*AZMAX1(Root2ndStructChemElmnt_pvr(NE,imycorrhz,LL,NR,NZ))*FWODRE(NE,k_woody_litr)

                        LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,LL,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,LL,NZ)+CFOPE(NE,iroot,M,NZ) &
                          *FSNCM*AZMAX1(Root2ndStructChemElmnt_pvr(NE,imycorrhz,LL,NR,NZ))*FWODRE(NE,k_fine_litr)

                        LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,LL,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,LL,NZ)+CFOPE(NE,inonstruct,M,NZ) &
                          *FSNCP*AZMAX1(RootMycoNonstructElmnt_vr(NE,imycorrhz,LL,NZ))
                      ENDDO D6450
                      Root2ndStructChemElmnt_pvr(NE,imycorrhz,LL,NR,NZ)=AZMAX1(Root2ndStructChemElmnt_pvr(NE,imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)

                      RootMycoNonstructElmnt_vr(NE,imycorrhz,LL,NZ)=AZMAX1(RootMycoNonstructElmnt_vr(NE,imycorrhz,LL,NZ))*(1.0_r8-FSNCP)

                    ENDDO
                    SecndRootLen(imycorrhz,LL,NR,NZ)=AZMAX1(SecndRootLen(imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)
                  ENDIF
                ENDDO D5105
              ENDIF
!
              call PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,RootCYieldO2ltd,RootNLigthSatCarboxyRate_nodewthElmnt,&
                Root1stExtension,TotPrimRootLen,Root1stC)
            ENDIF
!
!
            IF(L.EQ.NIXBotRootLayer_rpft(NR,NZ))THEN
              call WithdrawPrimRoot(L,NR,NZ,N,RTDP1X,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootAreaPopu)
            ENDIF

!
!     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
!
            IF(Root1stStructChemElmnt_pvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
               RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)+Root1stStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
              Root1stStructChemElmnt_pvr(ielmc,N,L,NR,NZ)=0._r8
            ENDIF
            IF(Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
               RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)+Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
              Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)=0._r8
            ENDIF
!
!     TOTAL PRIMARY ROOT LENGTH AND MASS
!
!     TotPrimRootLen=total primary root length
!     Root1stC=total primary root C mass
!     PrimRootLen=primary root length in soil layer
!     WTRT1=primary root C mass in soil layer
!     NIXBotRootLayer_rpft=deepest root layer
!
            TotPrimRootLen=TotPrimRootLen+PrimRootLen(N,L,NR,NZ)
            Root1stC=Root1stC+Root1stStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
            NIXBotRootLayer_rpft(NR,NZ)=MIN(NIXBotRootLayer_rpft(NR,NZ),MaxNumRootLays)
            IF(L.EQ.NIXBotRootLayer_rpft(NR,NZ))NRX(N,NR)=1
          ENDIF
        ENDIF
        TotPrimRootLen=TotPrimRootLen+PrimRootLen(N,L,NR,NZ)
        Root1stC=Root1stC+Root1stStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
!     ENDIF
      ENDIF
      NIXBotRootLayer_pft(NZ)=MAX(NIXBotRootLayer_pft(NZ),NIXBotRootLayer_rpft(NR,NZ))
  ENDDO D5050
  end associate
  end subroutine GrowRootAxes

!------------------------------------------------------------------------------------------

  subroutine PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

  implicit none
  integer, intent(in) :: N,L,NZ,NR
  real(r8), intent(in) :: RCO2X,RCO2XM
  real(r8), intent(out) :: FSNC1
  integer :: M,NE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: RootChemElmntRemob(NumOfPlantChemElmnts)
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SNCR,SNCRM
! begin_execution
  associate(                          &
    Root1stChemElmnt  =>  plt_biom%Root1stChemElmnt     , &
    RootNonstructElementConcpft_vr  =>  plt_biom%RootNonstructElementConcpft_vr     , &
    ZEROP   =>  plt_biom%ZEROP      , &
    ZERO    =>  plt_site%ZERO       , &
    icwood  =>  pltpar%icwood       , &
    k_woody_litr=> pltpar%k_woody_litr,&
    k_fine_litr=> pltpar%k_fine_litr, &
    iroot   =>  pltpar%iroot        , &
    FWODRE  =>  plt_allom%FWODRE    , &
    CFOPE   =>  plt_soilchem%CFOPE  , &
    RootAutoRO2Limiter_pvr    =>  plt_rbgc%RootAutoRO2Limiter_pvr       , &
    LitterFallChemElmnt_pftvr    =>  plt_bgcr%LitterFallChemElmnt_pftvr       , &
    iPlantCalendar_brch  =>  plt_pheno%iPlantCalendar_brch    , &
    NumOfMainBranch_pft     =>  plt_morph%NumOfMainBranch_pft         &
  )
!
!     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
  IF(iPlantCalendar_brch(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).NE.0 &
    .AND.RootNonstructElementConcpft_vr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElementConcpft_vr(ielmn,N,L,NZ),RootNonstructElementConcpft_vr(ielmn,N,L,NZ)+RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CNKI) &
      ,safe_adb(RootNonstructElementConcpft_vr(ielmp,N,L,NZ),RootNonstructElementConcpft_vr(ielmp,N,L,NZ)+RootNonstructElementConcpft_vr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElementConcpft_vr(ielmc,N,L,NZ),RootNonstructElementConcpft_vr(ielmc,N,L,NZ)+RootNonstructElementConcpft_vr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElementConcpft_vr(ielmc,N,L,NZ),RootNonstructElementConcpft_vr(ielmc,N,L,NZ)+RootNonstructElementConcpft_vr(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC=0._r8
    CNC=0._r8
    CPC=0._r8
  ENDIF
  RCCC=RCCZR+CCC*RCCYR
  RCCN=CNC*RCCXR
  RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P DURING PRIMARY ROOT REMOBILIZATION
!     DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RootAutoRO2Limiter_pvr=constraint by O2 consumption on all root processes
!     RootChemElmntRemob(ielmc),RootChemElmntRemob(ielmn),RootChemElmntRemob(ielmp)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC1=fraction of primary root C to be remobilized
!
  IF(-RCO2XM.GT.0.0_r8)THEN
    IF(-RCO2XM.LT.Root1stChemElmnt(ielmc,N,NR,NZ)*RCCC)THEN
      SNCRM=-RCO2XM
    ELSE
      SNCRM=AZMAX1(Root1stChemElmnt(ielmc,N,NR,NZ)*RCCC)
    ENDIF
  ELSE
    SNCRM=0._r8
  ENDIF
  IF(-RCO2X.GT.0.0_r8)THEN
    IF(-RCO2X.LT.Root1stChemElmnt(ielmc,N,NR,NZ)*RCCC)THEN
      SNCR=-RCO2X
    ELSE
      SNCR=AZMAX1(Root1stChemElmnt(ielmc,N,NR,NZ)*RCCC)*RootAutoRO2Limiter_pvr(N,L,NZ)
    ENDIF
  ELSE
    SNCR=0._r8
  ENDIF
  IF(SNCR.GT.0.0_r8.AND.Root1stChemElmnt(ielmc,N,NR,NZ).GT.ZEROP(NZ))THEN
    RootChemElmntRemob(ielmc)=RCCC*Root1stChemElmnt(ielmc,N,NR,NZ)
    RootChemElmntRemob(ielmn)=Root1stChemElmnt(ielmn,N,NR,NZ)*(RCCN+(1.0_r8-RCCN) &
      *RootChemElmntRemob(ielmc)/Root1stChemElmnt(ielmc,N,NR,NZ))
    RootChemElmntRemob(ielmp)=Root1stChemElmnt(ielmp,N,NR,NZ)*(RCCP+(1.0_r8-RCCP) &
      *RootChemElmntRemob(ielmc)/Root1stChemElmnt(ielmc,N,NR,NZ))
    IF(RootChemElmntRemob(ielmc).GT.ZEROP(NZ))THEN
      FSNC1=AZMAX1(AMIN1(1.0_r8,SNCR/RootChemElmntRemob(ielmc)))
    ELSE
      FSNC1=1.0_r8
    ENDIF
  ELSE
    RootChemElmntRemob(1:NumOfPlantChemElmnts)=0._r8
    FSNC1=0._r8
  ENDIF
!
!     PRIMARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RootChemElmntRemob(ielmc),RootChemElmntRemob(ielmn),RootChemElmntRemob(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
  D6355: DO M=1,jsken
    DO NE=1,NumOfPlantChemElmnts    
      LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ) &
        +CFOPE(NE,icwood,M,NZ)*FSNC1*(Root1stChemElmnt(NE,N,NR,NZ)-RootChemElmntRemob(NE))*FWODRE(NE,k_woody_litr)

      LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ) &
        +CFOPE(NE,iroot,M,NZ)*FSNC1*(Root1stChemElmnt(NE,N,NR,NZ)-RootChemElmntRemob(NE))*FWODRE(NE,k_fine_litr)
    ENDDO
  ENDDO D6355

  end associate
  end subroutine PrimRootRemobilization

!------------------------------------------------------------------------------------------
  subroutine PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,RootCYieldO2ltd,RootNLigthSatCarboxyRate_nodewthElmnt,&
    Root1stExtension,TotPrimRootLen,Root1stC)
  implicit none
  integer, intent(in) :: L,L1,N,NR,NZ
  real(r8), intent(in):: WFNR,FRTN,RootCYieldO2ltd
  real(r8), intent(in) :: RootNLigthSatCarboxyRate_nodewthElmnt(NumOfPlantChemElmnts)
  real(r8), intent(inout) :: TotPrimRootLen,Root1stC
  real(r8), intent(out):: Root1stExtension
  real(r8) :: FGROL,FGROZ
  integer :: NE
  REAL(R8) :: XFRE(NumOfPlantChemElmnts)
! begin_execution
  associate(                             &
    Root1stStructChemElmnt_pvr    =>  plt_biom%Root1stStructChemElmnt_pvr      , &
    Root1stChemElmnt    =>  plt_biom%Root1stChemElmnt      , &
    RootProteinC_pvr     =>  plt_biom%RootProteinC_pvr       , &
     PopuPlantRootC_vr    =>  plt_biom% PopuPlantRootC_vr      , &
     RootMycoNonstructElmnt_vr    =>  plt_biom%RootMycoNonstructElmnt_vr      , &
    ZEROP     =>  plt_biom%ZEROP       , &
    rCNNonstructRemob_pft     =>  plt_allom%rCNNonstructRemob_pft      , &
    rCPNonstructRemob_pft      =>  plt_allom%rCPNonstructRemob_pft       , &
    FWODRE    =>  plt_allom%FWODRE     , &
    CumSoilThickness    =>  plt_site%CumSoilThickness      , &
    DLYR3     =>  plt_site%DLYR3       , &
    MaxNumRootLays        =>  plt_site%MaxNumRootLays          , &
    PlantPopulation_pft        =>  plt_site%PlantPopulation_pft          , &
    PSIRootTurg_vr     =>  plt_ew%PSIRootTurg_vr         , &
    PSIRootOSMO_vr     =>  plt_ew%PSIRootOSMO_vr         , &
    PSIRoot_vr    =>  plt_ew%PSIRoot_vr        , &
    k_woody_litr=> pltpar%k_woody_litr , &
    k_fine_litr=> pltpar%k_fine_litr   , &
    PrimRootSpecLen    =>  plt_morph%PrimRootSpecLen     , &
    PrimRootRadius    =>  plt_morph%PrimRootRadius     , &
    PrimRootLen     =>  plt_morph%PrimRootLen      , &
    PrimRootDepth    =>  plt_morph%PrimRootDepth     , &
    NGTopRootLayer_pft       =>  plt_morph%NGTopRootLayer_pft        , &
    NIXBotRootLayer_rpft      =>  plt_morph%NIXBotRootLayer_rpft       , &
    SeedinDepth     =>  plt_morph%SeedinDepth        &
  )
!     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     Root1stExtension=primary root length extension
!     RootCYieldO2ltd=primary root C growth ltd by O2
!     PrimRootSpecLen=specific primary root length from startq.f
!     PP=PFT population
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     RootNLigthSatCarboxyRate_nodewthElmnt(ielmc),RootNLigthSatCarboxyRate_nodewthElmnt(ielmn),RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=net primary root C,N,P growth
!     RTDP1=primary root depth from soil surface
!     SeedinDepth=seeding depth
!     FSNC1=fraction of primary root C to be remobilized
!     PrimRootLen=primary root length
!     RootNLigthSatCarboxyRate_nodewthElmnt(ielmc),RootNLigthSatCarboxyRate_nodewthElmnt(ielmn),RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=net root C,N,P growth
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     DLYR=soil layer thickness
!
  IF(RootNLigthSatCarboxyRate_nodewthElmnt(ielmc).LT.0.0.AND.Root1stChemElmnt(ielmc,N,NR,NZ).GT.ZEROP(NZ))THEN
    Root1stExtension=RootCYieldO2ltd*PrimRootSpecLen(N,NZ)/PlantPopulation_pft(NZ)*WFNR*FWODRE(ielmc,k_fine_litr) &
      +RootNLigthSatCarboxyRate_nodewthElmnt(ielmc)*(PrimRootDepth(N,NR,NZ)-SeedinDepth(NZ))/Root1stChemElmnt(ielmc,N,NR,NZ)
  ELSE
    Root1stExtension=RootCYieldO2ltd*PrimRootSpecLen(N,NZ)/PlantPopulation_pft(NZ)*WFNR*FWODRE(ielmc,k_fine_litr)
  ENDIF
  IF(L.LT.MaxNumRootLays)THEN
    Root1stExtension=AMIN1(DLYR3(L1),Root1stExtension)
  ENDIF
!
!     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
!     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
!     BOUNDARY OF CURRENT LAYER
!
!     Root1stExtension=primary root length extension
!     FGROL,FGROZ=fraction of Root1stExtension in current,next lower soil layer
!
  IF(Root1stExtension.GT.ZEROP(NZ).AND.L.LT.MaxNumRootLays)THEN
    FGROL=AZMAX1(AMIN1(1.0_r8,(CumSoilThickness(L)-PrimRootDepth(N,NR,NZ))/Root1stExtension))
    IF(FGROL.LT.1.0_r8)FGROL=0._r8
    FGROZ=AZMAX1(1.0_r8-FGROL)
  ELSE
    FGROL=1.0_r8
    FGROZ=0._r8
  ENDIF
!
!     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
!     AND AXIS NUMBER
!
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RootNLigthSatCarboxyRate_nodewthElmnt(ielmc),RootNLigthSatCarboxyRate_nodewthElmnt(ielmn),RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=net root C,N,P growth
!     Root1stExtension=primary root length extension
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     FGROL,FGROZ=fraction of Root1stExtension in current,next lower soil layer
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!     PrimRootLen=primary root length
!
  PrimRootDepth(N,NR,NZ)=PrimRootDepth(N,NR,NZ)+Root1stExtension

  DO NE=1,NumOfPlantChemElmnts
    Root1stChemElmnt(NE,N,NR,NZ)=Root1stChemElmnt(NE,N,NR,NZ)+RootNLigthSatCarboxyRate_nodewthElmnt(NE)
    Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ)=Root1stStructChemElmnt_pvr(NE,N,L,NR,NZ)+RootNLigthSatCarboxyRate_nodewthElmnt(NE)*FGROL
  ENDDO
  RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)+AMIN1(rCNNonstructRemob_pft(NZ)*Root1stStructChemElmnt_pvr(ielmn,N,L,NR,NZ) &
    ,rCPNonstructRemob_pft(NZ)*Root1stStructChemElmnt_pvr(ielmp,N,L,NR,NZ))
  PrimRootLen(N,L,NR,NZ)=PrimRootLen(N,L,NR,NZ)+Root1stExtension*FGROL
!
!     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
!     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
!     OF CURRENT SOIL LAYER
!
!     FGROZ=fraction of Root1stExtension in next lower soil layer
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     RootNLigthSatCarboxyRate_nodewthElmnt(ielmc),RootNLigthSatCarboxyRate_nodewthElmnt(ielmn),RootNLigthSatCarboxyRate_nodewthElmnt(ielmp)=net root C,N,P growth
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!     WTRTD=root C mass
!     PrimRootLen=primary root length
!     Root1stExtension=primary root length extension
!     FRTN=fraction of primary root sink strength in axis
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     PSIRoot_vr,PSIRootTurg_vr,PSIRootOSMO_vr=root total,turgor,osmotic water potential
!     NIXBotRootLayer_rpft=deepest root layer
!
  IF(FGROZ.GT.0.0_r8)THEN
    DO NE=1,NumOfPlantChemElmnts
      Root1stStructChemElmnt_pvr(NE,N,L1,NR,NZ)=Root1stStructChemElmnt_pvr(NE,N,L1,NR,NZ)+RootNLigthSatCarboxyRate_nodewthElmnt(NE)*FGROZ
    ENDDO
    RootProteinC_pvr(N,L1,NZ)=RootProteinC_pvr(N,L1,NZ)+AMIN1(rCNNonstructRemob_pft(NZ)*Root1stStructChemElmnt_pvr(ielmn,N,L1,NR,NZ) &
      ,rCPNonstructRemob_pft(NZ)*Root1stStructChemElmnt_pvr(ielmp,N,L1,NR,NZ))
    PopuPlantRootC_vr(N,L1,NZ)= PopuPlantRootC_vr(N,L1,NZ)+Root1stStructChemElmnt_pvr(ielmc,N,L1,NR,NZ)
    PrimRootLen(N,L1,NR,NZ)=PrimRootLen(N,L1,NR,NZ)+Root1stExtension*FGROZ
    PrimRootRadius(N,L1,NZ)=PrimRootRadius(N,L,NZ)
    TotPrimRootLen=TotPrimRootLen+PrimRootLen(N,L1,NR,NZ)
    Root1stC=Root1stC+Root1stStructChemElmnt_pvr(ielmc,N,L1,NR,NZ)

    DO NE=1,NumOfPlantChemElmnts
      XFRE(NE)=FRTN* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
      RootMycoNonstructElmnt_vr(NE,N,L,NZ)= RootMycoNonstructElmnt_vr(NE,N,L,NZ)-XFRE(NE)
      RootMycoNonstructElmnt_vr(NE,N,L1,NZ)= RootMycoNonstructElmnt_vr(NE,N,L1,NZ)+XFRE(NE)
    ENDDO
    PSIRoot_vr(N,L1,NZ)=PSIRoot_vr(N,L,NZ)
    PSIRootOSMO_vr(N,L1,NZ)=PSIRootOSMO_vr(N,L,NZ)
    PSIRootTurg_vr(N,L1,NZ)=PSIRootTurg_vr(N,L,NZ)
    NIXBotRootLayer_rpft(NR,NZ)=MAX(NGTopRootLayer_pft(NZ),L+1)
  ENDIF
  end associate
  end subroutine PrimRootExtension
!------------------------------------------------------------------------------------------

  subroutine WithdrawPrimRoot(L,NR,NZ,N,RTDP1X,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootAreaPopu)
  implicit none
  integer, intent(in) :: L,NR,NZ,N
  real(r8), intent(in):: RTDP1X
  real(r8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,10),Root2ndSink_pvr(jroots,JZ1,10)
  real(r8), intent(in) :: RootAreaPopu
  integer :: LL,NN,NE,NTG
  real(r8) :: XFRD,XFRW,FRTN
  real(r8) :: XFRE(NumOfPlantChemElmnts)

! begin_execution
  associate(                             &
    Root1stStructChemElmnt_pvr   =>  plt_biom%Root1stStructChemElmnt_pvr   , &
    Root2ndStructChemElmnt_pvr   =>  plt_biom%Root2ndStructChemElmnt_pvr   , &
    RootMycoNonstructElmnt_vr   =>  plt_biom%RootMycoNonstructElmnt_vr   , &
    RootProteinC_pvr    =>  plt_biom%RootProteinC_pvr    , &
    PopuPlantRootC_vr   =>  plt_biom% PopuPlantRootC_vr   , &
    RootNodueChemElmnt_pvr   =>  plt_biom%RootNodueChemElmnt_pvr   , &
    ZEROP    =>  plt_biom%ZEROP    , &
    RootNoduleNonstructElmnt_vr  =>  plt_biom%RootNoduleNonstructElmnt_vr  , &
    RootGasLoss_disturb    =>  plt_bgcr%RootGasLoss_disturb    , &
    CumSoilThickness   =>  plt_site%CumSoilThickness   , &
    ZEROS2   =>  plt_site%ZEROS2   , &
    DLYR3    =>  plt_site%DLYR3    , &
    trcg_rootml     =>  plt_rbgc%trcg_rootml , &
    trcs_rootml => plt_rbgc%trcs_rootml, &
    VLSoilPoreMicP     =>  plt_soilchem%VLSoilPoreMicP , &
    PrimRootXNumL     =>  plt_morph%PrimRootXNumL    , &
    SecndRootXNum_pvr     =>  plt_morph%SecndRootXNum_pvr    , &
    PrimRootLen    =>  plt_morph%PrimRootLen   , &
    SecndRootXNum_rpvr    =>  plt_morph%SecndRootXNum_rpvr   , &
    PrimRootDepth   =>  plt_morph%PrimRootDepth  , &
    NGTopRootLayer_pft      =>  plt_morph%NGTopRootLayer_pft     , &
    iPlantNfixType   =>  plt_morph%iPlantNfixType  , &
    MY       =>  plt_morph%MY      , &
    SeedinDepth    =>  plt_morph%SeedinDepth   , &
    NIXBotRootLayer_rpft     =>  plt_morph%NIXBotRootLayer_rpft      &
  )
!     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
!     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
!     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
!     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
!
!     NIXBotRootLayer_rpft=deepest root layer
!     VLSoilPoreMicP=soil layer volume excluding macropore, rocks
!     RTDP1X=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     SeedinDepth=seeding depth
!     FRTN=fraction of primary root sink strength in axis
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     PrimRootLen=primary root length
!     RootProteinC_pvr=root protein C mass
!     WTRTD=root C mass
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root

  D5115: DO LL=L,NGTopRootLayer_pft(NZ)+1,-1
    IF(VLSoilPoreMicP(LL-1).GT.ZEROS2.AND.(RTDP1X.LT.CumSoilThickness(LL-1).OR.RTDP1X.LT.SeedinDepth(NZ)))THEN
      IF(RootSinkC_vr(N,LL).GT.ZEROP(NZ))THEN
        FRTN=(Root1stSink_pvr(N,LL,NR)+Root2ndSink_pvr(N,LL,NR))/RootSinkC_vr(N,LL)
      ELSE
        FRTN=1.0_r8
      ENDIF
      D5110: DO NN=1,MY(NZ)
        PrimRootLen(NN,LL-1,NR,NZ)=PrimRootLen(NN,LL-1,NR,NZ)+PrimRootLen(NN,LL,NR,NZ)
        PrimRootLen(NN,LL,NR,NZ)=0._r8
        DO NE=1,NumOfPlantChemElmnts
          Root1stStructChemElmnt_pvr(NE,NN,LL-1,NR,NZ)=Root1stStructChemElmnt_pvr(NE,NN,LL-1,NR,NZ)+Root1stStructChemElmnt_pvr(NE,NN,LL,NR,NZ)

          Root2ndStructChemElmnt_pvr(NE,NN,LL-1,NR,NZ)=Root2ndStructChemElmnt_pvr(NE,NN,LL-1,NR,NZ)+Root2ndStructChemElmnt_pvr(NE,NN,LL,NR,NZ)

          Root1stStructChemElmnt_pvr(NE,NN,LL,NR,NZ)=0._r8

          Root2ndStructChemElmnt_pvr(NE,NN,LL,NR,NZ)=0._r8

          XFRE(NE)=FRTN* RootMycoNonstructElmnt_vr(NE,NN,LL,NZ)

           RootMycoNonstructElmnt_vr(NE,NN,LL,NZ)= RootMycoNonstructElmnt_vr(NE,NN,LL,NZ)-XFRE(NE)

           RootMycoNonstructElmnt_vr(NE,NN,LL-1,NZ)= RootMycoNonstructElmnt_vr(NE,NN,LL-1,NZ)+XFRE(NE)
        ENDDO
        XFRW=FRTN*RootProteinC_pvr(NN,L,NZ)
        XFRD=FRTN* PopuPlantRootC_vr(NN,LL,NZ)

        RootProteinC_pvr(NN,LL,NZ)=RootProteinC_pvr(NN,LL,NZ)-XFRW
         PopuPlantRootC_vr(NN,LL,NZ)= PopuPlantRootC_vr(NN,LL,NZ)-XFRD
        RootProteinC_pvr(NN,LL-1,NZ)=RootProteinC_pvr(NN,LL-1,NZ)+XFRW
         PopuPlantRootC_vr(NN,LL-1,NZ)= PopuPlantRootC_vr(NN,LL-1,NZ)+XFRD
!
!     WITHDRAW GASES IN PRIMARY ROOTS
!
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
!     FRTN=fraction of primary root sink strength in axis
!
        DO NTG=idg_beg,idg_end-1
          RootGasLoss_disturb(NTG,NZ)=RootGasLoss_disturb(NTG,NZ)-FRTN &
            *(trcg_rootml(idg_CO2,NN,LL,NZ)+trcs_rootml(idg_CO2,NN,LL,NZ))
          trcg_rootml(NTG,NN,LL,NZ)=(1.0_r8-FRTN)*trcg_rootml(NTG,NN,LL,NZ)
          trcs_rootml(NTG,NN,LL,NZ)=(1.0_r8-FRTN)*trcs_rootml(NTG,NN,LL,NZ)
        ENDDO
      ENDDO D5110
!
!     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
!
!     SecndRootXNum_rpvr,SecndRootXNum_pvr=number of secondary root axes
!     PrimRootXNumL=number of primary root axes
!     PrimRootLen=primary root length
!     CumSoilThickness=depth from soil surface to layer bottom
!     SeedinDepth=seeding depth
!
      SecndRootXNum_pvr(N,LL,NZ)=SecndRootXNum_pvr(N,LL,NZ)-SecndRootXNum_rpvr(N,LL,NR,NZ)
      SecndRootXNum_pvr(N,LL-1,NZ)=SecndRootXNum_pvr(N,LL-1,NZ)+SecndRootXNum_rpvr(N,LL,NR,NZ)
      SecndRootXNum_rpvr(N,LL,NR,NZ)=0._r8
      PrimRootXNumL(N,LL,NZ)=PrimRootXNumL(N,LL,NZ)-RootAreaPopu
      IF(LL-1.GT.NGTopRootLayer_pft(NZ))THEN
        PrimRootLen(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness(LL-1)-PrimRootDepth(N,NR,NZ))
      ELSE
        PrimRootLen(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness(LL-1)-PrimRootDepth(N,NR,NZ)) &
          -(SeedinDepth(NZ)-CumSoilThickness(LL-2))
      ENDIF
!
!     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
!
!     iPlantNfixType=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     FRTN=fraction of primary root sink strength in axis
!     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
!
      IF(is_root_N2fix(iPlantNfixType(NZ)))THEN
        DO NE=1,NumOfPlantChemElmnts
          XFRE(NE)=FRTN*RootNodueChemElmnt_pvr(NE,LL,NZ)
          RootNodueChemElmnt_pvr(NE,LL,NZ)=RootNodueChemElmnt_pvr(NE,LL,NZ)-XFRE(NE)
          RootNodueChemElmnt_pvr(NE,LL-1,NZ)=RootNodueChemElmnt_pvr(NE,LL-1,NZ)+XFRE(NE)
          XFRE(NE)=FRTN*RootNoduleNonstructElmnt_vr(NE,LL,NZ)
          RootNoduleNonstructElmnt_vr(NE,LL,NZ)=RootNoduleNonstructElmnt_vr(NE,LL,NZ)-XFRE(NE)
          RootNoduleNonstructElmnt_vr(NE,LL-1,NZ)=RootNoduleNonstructElmnt_vr(NE,LL-1,NZ)+XFRE(NE)
        ENDDO
      ENDIF
      NIXBotRootLayer_rpft(NR,NZ)=MAX(NGTopRootLayer_pft(NZ),LL-1)
    ELSE
      EXIT
    ENDIF
  ENDDO D5115
  end associate
  end subroutine WithdrawPrimRoot
!------------------------------------------------------------------------------------------

  subroutine NonstructlBiomTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,BegRemoblize)
  implicit none
  integer, intent(in) :: I,J,NZ,BegRemoblize
  real(r8), intent(in):: PTRT
  real(r8), INTENT(IN) :: RootSinkC_vr(2,JZ1)
  real(r8),intent(in) :: Root1stSink_pvr(2,JZ1,10),Root2ndSink_pvr(2,JZ1,10)
  real(r8),intent(in) :: RootSinkC(2)
  integer :: L,NB,N,NR,NE
  real(r8) :: ZPOOLS,ZPOOLT
  real(r8) :: ZPOOLB
  real(r8) :: ZPOOLD,EPOOLD
  real(r8) :: XFRPX,XFRCX,XFRNX
  real(r8) :: WTLSBZ(NumOfCanopyLayers1)
  real(r8) :: CPOOLZ(NumOfCanopyLayers1),ZPOOLZ(NumOfCanopyLayers1),PPOOLZ(NumOfCanopyLayers1)
  REAL(R8) :: FWTR(JZ1),FWTB(JP1)
  real(r8) :: CPOOLT
  real(r8) :: CNL,CPL,CPOOLD
  real(r8) :: CPOOLB,CPOOLS
  real(r8) :: FWTC
  real(r8) :: FWTS
  real(r8) :: PPOOLB
  real(r8) :: PPOOLD
  real(r8) :: PPOOLT
  real(r8) :: PTSHTR
  real(r8) :: PPOOLS
  real(r8) :: WTPLTT
  real(r8) :: WTRTLX
  real(r8) :: WTRTD1
  real(r8) :: WTSTKT
  real(r8) :: WTRSVT,WTRSNT,WTRSPT
  real(r8) :: WTRSVD,WTRSND,WTRSPD
  real(r8) :: WTRTD2,WTLSBX,WTLSBB
  real(r8) :: WTRTLR
  real(r8) :: XFRE(NumOfPlantChemElmnts)
!     begin_execution
  associate(                                &
    FWODBE     =>   plt_allom%FWODBE  , &
    FWODRE     =>   plt_allom%FWODRE  , &
    RootNonstructElementConcpft_vr     =>   plt_biom%RootNonstructElementConcpft_vr   , &
    Root1stChemElmnt     =>   plt_biom%Root1stChemElmnt   , &
    Root1stStructChemElmnt_pvr     =>   plt_biom%Root1stStructChemElmnt_pvr   , &
    Root2ndStructChemElmnt_pvr     =>   plt_biom%Root2ndStructChemElmnt_pvr   , &
    RootStructBiomC_vr     =>   plt_biom%RootStructBiomC_vr   , &
     RootMycoNonstructElmnt_vr     =>   plt_biom%RootMycoNonstructElmnt_vr   , &
    LeafPetioleBiomassC_brch    =>   plt_biom%LeafPetioleBiomassC_brch   , &
    NonstructElmnt_brch     =>   plt_biom%NonstructElmnt_brch   , &
     PopuPlantRootC_vr     =>   plt_biom% PopuPlantRootC_vr   , &
    StalkBiomassC_brch     =>   plt_biom%StalkBiomassC_brch   , &
    ReserveChemElmnts_brch    =>   plt_biom%ReserveChemElmnts_brch  , &
    NonstructalChemElmnts_pft      =>   plt_biom%NonstructalChemElmnts_pft    , &
    CanopyLeafShethC_pft       =>   plt_biom%CanopyLeafShethC_pft     , &
    RootChemElmnts_pft      =>   plt_biom%RootChemElmnts_pft    , &
    ZEROL      =>   plt_biom%ZEROL    , &
    ZEROP      =>   plt_biom%ZEROP    , &
    iPlantTurnoverPattern     =>   plt_pheno%iPlantTurnoverPattern  , &
    iPlantMorphologyType     =>   plt_pheno%iPlantMorphologyType  , &
    iPlantBranchState      =>   plt_pheno%iPlantBranchState   , &
    iPlantPhenologyPattern_pft     =>   plt_pheno%iPlantPhenologyPattern_pft  , &
    ShutRutNonstructElmntConducts     =>   plt_pheno%ShutRutNonstructElmntConducts  , &
    iPlantCalendar_brch     =>   plt_pheno%iPlantCalendar_brch  , &
    HourCounter4LeafOut_brch       =>   plt_pheno%HourCounter4LeafOut_brch    , &
    RCO2A      =>   plt_rbgc%RCO2A    , &
    GrossResp_pft      =>   plt_bgcr%GrossResp_pft    , &
    ECO_ER_col       =>   plt_bgcr%ECO_ER_col     , &
    Eco_AutoR_col       =>   plt_bgcr%Eco_AutoR_col     , &
    NU         =>   plt_site%NU       , &
    ZERO       =>   plt_site%ZERO     , &
    k_woody_litr=> pltpar%k_woody_litr,&
    k_fine_litr=> pltpar%k_fine_litr, &
    NIXBotRootLayer_pft       =>   plt_morph%NIXBotRootLayer_pft    , &
    NIXBotRootLayer_rpft       =>   plt_morph%NIXBotRootLayer_rpft    , &
    SecndRootRadius      =>   plt_morph%SecndRootRadius   , &
    NI         =>   plt_morph%NI      , &
    MY         =>   plt_morph%MY      , &
    NumRootAxes_pft       =>   plt_morph%NumRootAxes_pft    , &
    NumOfBranches_pft        =>   plt_morph%NumOfBranches_pft       &
  )
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH LEAVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!     WHEN SEASONAL STORAGE C IS NOT BEING MOBILIZED
!
!     iPlantBranchState=branch living flag: 0=alive,1=dead
!     HourCounter4LeafOut_brch=hourly leafout counter
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     LeafPetioleBiomassC_brch=leaf+petiole mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!
  IF(NumOfBranches_pft(NZ).GT.1)THEN
    WTPLTT=0._r8
    CPOOLT=0._r8
    ZPOOLT=0._r8
    PPOOLT=0._r8
    D300: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
        IF(HourCounter4LeafOut_brch(NB,NZ).GT.ATRPX(iPlantPhenologyPattern_pft(NZ)))THEN
          WTLSBZ(NB)=AZMAX1(LeafPetioleBiomassC_brch(NB,NZ))
          CPOOLZ(NB)=AZMAX1(NonstructElmnt_brch(ielmc,NB,NZ))
          ZPOOLZ(NB)=AZMAX1(NonstructElmnt_brch(ielmn,NB,NZ))
          PPOOLZ(NB)=AZMAX1(NonstructElmnt_brch(ielmp,NB,NZ))
          WTPLTT=WTPLTT+WTLSBZ(NB)
          CPOOLT=CPOOLT+CPOOLZ(NB)
          ZPOOLT=ZPOOLT+ZPOOLZ(NB)
          PPOOLT=PPOOLT+PPOOLZ(NB)
        ENDIF
      ENDIF
    ENDDO D300
    D305: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
        IF(HourCounter4LeafOut_brch(NB,NZ).GT.ATRPX(iPlantPhenologyPattern_pft(NZ)))THEN
          IF(WTPLTT.GT.ZEROP(NZ).AND.CPOOLT.GT.ZEROP(NZ))THEN
            CPOOLD=CPOOLT*WTLSBZ(NB)-CPOOLZ(NB)*WTPLTT
            ZPOOLD=ZPOOLT*CPOOLZ(NB)-ZPOOLZ(NB)*CPOOLT
            PPOOLD=PPOOLT*CPOOLZ(NB)-PPOOLZ(NB)*CPOOLT
            XFRE(ielmc)=0.01_r8*CPOOLD/WTPLTT
            XFRE(ielmn)=0.01_r8*ZPOOLD/CPOOLT
            XFRE(ielmp)=0.01_r8*PPOOLD/CPOOLT
            DO NE=1,NumOfPlantChemElmnts
              NonstructElmnt_brch(NE,NB,NZ)=NonstructElmnt_brch(NE,NB,NZ)+XFRE(NE)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO D305
  ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH STALK RESERVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     iPlantBranchState=branch living flag: 0=alive,1=dead
!     StalkBiomassC_brch=stalk sapwood mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
!
  IF(NumOfBranches_pft(NZ).GT.1)THEN
    WTSTKT=0._r8
    WTRSVT=0._r8
    WTRSNT=0._r8
    WTRSPT=0._r8
    D330: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
        IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
          WTSTKT=WTSTKT+StalkBiomassC_brch(NB,NZ)
          WTRSVT=WTRSVT+ReserveChemElmnts_brch(ielmc,NB,NZ)
          WTRSNT=WTRSNT+ReserveChemElmnts_brch(ielmn,NB,NZ)
          WTRSPT=WTRSPT+ReserveChemElmnts_brch(ielmp,NB,NZ)
        ENDIF
      ENDIF
    ENDDO D330
    IF(WTSTKT.GT.ZEROP(NZ).AND.WTRSVT.GT.ZEROP(NZ))THEN
      D335: DO NB=1,NumOfBranches_pft(NZ)
        IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
          IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
            WTRSVD=WTRSVT*StalkBiomassC_brch(NB,NZ)-ReserveChemElmnts_brch(ielmc,NB,NZ)*WTSTKT
            XFRE(ielmc)=0.1_r8*WTRSVD/WTSTKT
            ReserveChemElmnts_brch(ielmc,NB,NZ)=ReserveChemElmnts_brch(ielmc,NB,NZ)+XFRE(ielmc)
            WTRSND=WTRSNT*ReserveChemElmnts_brch(ielmc,NB,NZ)-ReserveChemElmnts_brch(ielmn,NB,NZ)*WTRSVT
            XFRE(ielmn)=0.1_r8*WTRSND/WTRSVT
            ReserveChemElmnts_brch(ielmn,NB,NZ)=ReserveChemElmnts_brch(ielmn,NB,NZ)+XFRE(ielmn)
            WTRSPD=WTRSPT*ReserveChemElmnts_brch(ielmc,NB,NZ)-ReserveChemElmnts_brch(ielmp,NB,NZ)*WTRSVT
            XFRE(ielmp)=0.1_r8*WTRSPD/WTRSVT
            ReserveChemElmnts_brch(ielmp,NB,NZ)=ReserveChemElmnts_brch(ielmp,NB,NZ)+XFRE(ielmp)
          ENDIF
        ENDIF
      ENDDO D335
    ENDIF
  ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BWTWEEN ROOT AND MYCORRHIZAE
!     IN EACH ROOTED SOIL LAYER FROM NON-STRUCTURAL C,N,P
!     CONCENTRATION DIFFERENCES
!
!     MY=mycorrhizal:1=no,2=yes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in 1:root,2:mycorrhizae
!     WTRTD=1:root,2:mycorrhizal C mass
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     FMYC=rate constant for root-mycorrhizal C,N,P exchange (h-1)
!
  IF(MY(NZ).EQ.imycorr_arbu)THEN
    D425: DO L=NU,NIXBotRootLayer_pft(NZ)
      IF(RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ).GT.ZEROP(NZ).AND. PopuPlantRootC_vr(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
!root
        WTRTD1= PopuPlantRootC_vr(ipltroot,L,NZ)
        WTRTD2=AMIN1( PopuPlantRootC_vr(ipltroot,L,NZ),AMAX1(FSNK &
          * PopuPlantRootC_vr(ipltroot,L,NZ), PopuPlantRootC_vr(imycorrhz,L,NZ)))
        WTPLTT=WTRTD1+WTRTD2
        IF(WTPLTT.GT.ZEROP(NZ))THEN
          CPOOLD=(RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)*WTRTD2- RootMycoNonstructElmnt_vr(ielmc,imycorrhz,L,NZ)*WTRTD1)/WTPLTT
          XFRE(ielmc)=FMYC*CPOOLD
           RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
           RootMycoNonstructElmnt_vr(ielmc,imycorrhz,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,imycorrhz,L,NZ)+XFRE(ielmc)
          CPOOLT= RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)+ RootMycoNonstructElmnt_vr(ielmc,imycorrhz,L,NZ)

          IF(CPOOLT.GT.ZEROP(NZ))THEN
            DO NE=2,NumOfPlantChemElmnts
              EPOOLD=(RootMycoNonstructElmnt_vr(NE,ipltroot,L,NZ)* RootMycoNonstructElmnt_vr(ielmc,imycorrhz,L,NZ) &
                - RootMycoNonstructElmnt_vr(NE,imycorrhz,L,NZ)* RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(NE)=FMYC*EPOOLD
               RootMycoNonstructElmnt_vr(NE,ipltroot,L,NZ)= RootMycoNonstructElmnt_vr(NE,ipltroot,L,NZ)-XFRE(NE)
               RootMycoNonstructElmnt_vr(NE,imycorrhz,L,NZ)= RootMycoNonstructElmnt_vr(NE,imycorrhz,L,NZ)+XFRE(NE)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO D425
  ENDIF
!
!     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
!     IN PERENNIALS
!
  IF(BegRemoblize.EQ.1.AND.iPlantPhenologyPattern_pft(NZ).NE.iplt_annual)THEN
    D5545: DO N=1,MY(NZ)
      D5550: DO L=NU,NI(NZ)
        IF(RootNonstructElementConcpft_vr(ielmc,N,L,NZ).GT.ZERO)THEN
          CNL=RootNonstructElementConcpft_vr(ielmc,N,L,NZ)/(RootNonstructElementConcpft_vr(ielmc,N,L,NZ) &
            +RootNonstructElementConcpft_vr(ielmn,N,L,NZ)/CNKI)
          CPL=RootNonstructElementConcpft_vr(ielmc,N,L,NZ)/(RootNonstructElementConcpft_vr(ielmc,N,L,NZ) &
            +RootNonstructElementConcpft_vr(ielmp,N,L,NZ)/CPKI)
        ELSE
          CNL=0._r8
          CPL=0._r8
        ENDIF
        XFRCX=FXFR(iPlantTurnoverPattern(NZ))*AZMAX1(RootMycoNonstructElmnt_vr(ielmc,N,L,NZ))
        XFRNX=FXFR(iPlantTurnoverPattern(NZ))*AZMAX1(RootMycoNonstructElmnt_vr(ielmn,N,L,NZ))*(1.0_r8+CNL)
        XFRPX=FXFR(iPlantTurnoverPattern(NZ))*AZMAX1(RootMycoNonstructElmnt_vr(ielmp,N,L,NZ))*(1.0_r8+CPL)
        XFRE(ielmc)=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
        XFRE(ielmn)=AMIN1(XFRNX,XFRE(ielmc)*CNMX,XFRPX*CNMX/CPMN*0.5_r8)
        XFRE(ielmp)=AMIN1(XFRPX,XFRE(ielmc)*CPMX,XFRNX*CPMX/CNMN*0.5_r8)
        DO NE=1,NumOfPlantChemElmnts
           RootMycoNonstructElmnt_vr(NE,N,L,NZ)= RootMycoNonstructElmnt_vr(NE,N,L,NZ)-XFRE(NE)
          NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)+XFRE(NE)
        ENDDO

      ENDDO D5550
    ENDDO D5545
  ENDIF
!
!     ROOT AND NODULE TOTALS
!
!     RootStructBiomC_vr,WTRTD=active,actual root C mass
!     WTRT1,WTRT2=primary,secondary root C mass in soil layer
!     GrossResp_pft=total PFT respiration
!     RCO2A=total root respiration
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!
  D5445: DO N=1,MY(NZ)
    D5450: DO L=NU,NI(NZ)
      RootStructBiomC_vr(N,L,NZ)=0._r8
       PopuPlantRootC_vr(N,L,NZ)=0._r8
      D5460: DO NR=1,NumRootAxes_pft(NZ)
        RootStructBiomC_vr(N,L,NZ)=RootStructBiomC_vr(N,L,NZ)+Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
         PopuPlantRootC_vr(N,L,NZ)= PopuPlantRootC_vr(N,L,NZ)+Root2ndStructChemElmnt_pvr(ielmc,N,L,NR,NZ)+Root1stStructChemElmnt_pvr(ielmc,N,L,NR,NZ)
      ENDDO D5460
      GrossResp_pft(NZ)=GrossResp_pft(NZ)+RCO2A(N,L,NZ)
      ECO_ER_col=ECO_ER_col+RCO2A(N,L,NZ)
      Eco_AutoR_col=Eco_AutoR_col+RCO2A(N,L,NZ)
    ENDDO D5450

    DO  NR=1,NumRootAxes_pft(NZ)
      RootStructBiomC_vr(N,NIXBotRootLayer_rpft(NR,NZ),NZ)=RootStructBiomC_vr(N,NIXBotRootLayer_rpft(NR,NZ),NZ)+Root1stChemElmnt(ielmc,N,NR,NZ)
    ENDDO
  ENDDO D5445
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
!
!     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
!     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
!
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     WTLS,WTRT=total PFT leaf+petiole,root C mass
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     RootSinkC_vr,RootSinkC=root layer,root system sink strength
!
!     IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_perennial)THEN
  IF(CanopyLeafShethC_pft(NZ).GT.ZEROP(NZ))THEN
    FWTC=AMIN1(1.0_r8,0.667_r8*RootChemElmnts_pft(ielmc,NZ)/CanopyLeafShethC_pft(NZ))
  ELSE
    FWTC=1.0_r8
  ENDIF
  IF(RootChemElmnts_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
    FWTS=AMIN1(1.0_r8,CanopyLeafShethC_pft(NZ)/(0.667_r8*RootChemElmnts_pft(ielmc,NZ)))
  ELSE
    FWTS=1.0_r8
  ENDIF
!     ELSE
!     FWTC=1.0_r8
!     FWTS=1.0_r8
!     ENDIF
  D290: DO L=NU,NI(NZ)
    IF(RootSinkC(1).GT.ZEROP(NZ))THEN
      FWTR(L)=AZMAX1(RootSinkC_vr(ipltroot,L)/RootSinkC(1))
    ELSE
      FWTR(L)=1.0_r8
    ENDIF
  ENDDO D290
!     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
!     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
!
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!
  CanopyLeafShethC_pft(NZ)=0._r8
  D309: DO NB=1,NumOfBranches_pft(NZ)
    CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetioleBiomassC_brch(NB,NZ)
  ENDDO D309
!
!     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
!     OF TOTAL SINK STRENGTH OF THE CANOPY
!
!     iPlantBranchState=branch living flag: 0=alive,1=dead
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date for setting final seed number
!     FWTB=branch sink weighting factor
!     ShutRutNonstructElmntConducts=rate constant for equilibrating shoot-root nonstructural C concn from PFT file
!     PTRT=allocation to leaf+petiole used to modify ShutRutNonstructElmntConductsin annuals
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     FWODB=C woody fraction in branch:0=woody,1=non-woody
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  D310: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
      IF(CanopyLeafShethC_pft(NZ).GT.ZEROP(NZ))THEN
        FWTB(NB)=AZMAX1(LeafPetioleBiomassC_brch(NB,NZ)/CanopyLeafShethC_pft(NZ))
      ELSE
        FWTB(NB)=1.0_r8
      ENDIF
      IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual)THEN
        PTSHTR=ShutRutNonstructElmntConducts(NZ)*PTRT**0.167_r8
      ELSE
        PTSHTR=ShutRutNonstructElmntConducts(NZ)
      ENDIF
      D415: DO L=NU,NI(NZ)
        WTLSBX=LeafPetioleBiomassC_brch(NB,NZ)*FWODBE(ielmc,k_fine_litr)*FWTR(L)*FWTC
        WTRTLX=RootStructBiomC_vr(ipltroot,L,NZ)*FWODRE(ielmc,k_fine_litr)*FWTB(NB)*FWTS
        WTLSBB=AZMAX1(WTLSBX,FSNK*WTRTLX)
        WTRTLR=AZMAX1(WTRTLX,FSNK*WTLSBX)
        WTPLTT=WTLSBB+WTRTLR
        IF(WTPLTT.GT.ZEROP(NZ))THEN
          CPOOLB=AZMAX1(NonstructElmnt_brch(ielmc,NB,NZ)*FWTR(L))
          CPOOLS=AZMAX1(RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)*FWTB(NB))
          CPOOLD=(CPOOLB*WTRTLR-CPOOLS*WTLSBB)/WTPLTT
          XFRE(ielmc)=PTSHTR*CPOOLD
          CPOOLT=CPOOLS+CPOOLB
          IF(CPOOLT.GT.ZEROP(NZ))THEN
            ZPOOLB=AZMAX1(NonstructElmnt_brch(ielmn,NB,NZ)*FWTR(L))
            ZPOOLS=AZMAX1(RootMycoNonstructElmnt_vr(ielmn,ipltroot,L,NZ)*FWTB(NB))
            ZPOOLD=(ZPOOLB*CPOOLS-ZPOOLS*CPOOLB)/CPOOLT
            XFRE(ielmn)=PTSHTR*ZPOOLD
            PPOOLB=AZMAX1(NonstructElmnt_brch(ielmp,NB,NZ)*FWTR(L))
            PPOOLS=AZMAX1(RootMycoNonstructElmnt_vr(ielmp,ipltroot,L,NZ)*FWTB(NB))
            PPOOLD=(PPOOLB*CPOOLS-PPOOLS*CPOOLB)/CPOOLT
            XFRE(ielmp)=PTSHTR*PPOOLD
          ELSE
            XFRE(ielmn)=0._r8
            XFRE(ielmp)=0._r8
          ENDIF
          DO NE=1,NumOfPlantChemElmnts
            NonstructElmnt_brch(NE,NB,NZ)=NonstructElmnt_brch(NE,NB,NZ)-XFRE(NE)
             RootMycoNonstructElmnt_vr(NE,ipltroot,L,NZ)= RootMycoNonstructElmnt_vr(NE,ipltroot,L,NZ)+XFRE(NE)
          ENDDO

        ENDIF
      ENDDO D415
    ENDIF
  ENDDO D310
  end associate
  end subroutine NonstructlBiomTransfer

!------------------------------------------------------------------------------------------
  subroutine SummarizeRootSink(NZ,RootAreaPopu,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: RootAreaPopu
  real(r8),INTENT(OUT) :: RootSinkC_vr(jroots,JZ1)
  real(r8),intent(out) :: Root1stSink_pvr(jroots,JZ1,NumOfCanopyLayers1)
  real(r8),intent(out) :: Root2ndSink_pvr(jroots,JZ1,NumOfCanopyLayers1)
  real(r8),INTENT(OUT) :: RootSinkC(jroots)
  integer :: N,L,K,NR,NE
  REAL(R8) :: Root1stDepz_vr(NumOfCanopyLayers1,JZ1)
  real(r8) :: CUPRL,CUPRO,CUPRC
  real(r8) :: RTDPP,RTDPS,RTSKP
  real(r8) :: RTSKS

  associate(                              &
    RootMycoNonstructElmnt_vr     =>   plt_biom%RootMycoNonstructElmnt_vr    , &
    ZEROP      =>   plt_biom%ZEROP     , &
    iPlantMorphologyType     =>   plt_pheno%iPlantMorphologyType   , &
    VLSoilPoreMicP       =>   plt_soilchem%VLSoilPoreMicP  , &
    ZEROS2     =>   plt_site%ZEROS2    , &
    NU         =>   plt_site%NU        , &
    CumSoilThickness     =>   plt_site%CumSoilThickness    , &
    ZERO       =>   plt_site%ZERO      , &
    DLYR3      =>   plt_site%DLYR3     , &
    RUPH1B     =>   plt_rbgc%RUPH1B    , &
    RootH2PO4Uptake_pvr     =>   plt_rbgc%RootH2PO4Uptake_pvr    , &
    RUPH2B     =>   plt_rbgc%RUPH2B    , &
    RUOH2B     =>   plt_rbgc%RUOH2B    , &
    RUOH1P     =>   plt_rbgc%RUOH1P    , &
    RUCH1B     =>   plt_rbgc%RUCH1B    , &
    RUCH2B     =>   plt_rbgc%RUCH2B    , &
    RUONH4     =>   plt_rbgc%RUONH4    , &
    RUCH1P     =>   plt_rbgc%RUCH1P    , &
    RUCH2P     =>   plt_rbgc%RUCH2P    , &
    RUCNOB     =>   plt_rbgc%RUCNOB    , &
    RUCNO3     =>   plt_rbgc%RUCNO3    , &
    RDFOME     =>   plt_rbgc%RDFOME    , &
    RUPNOB     =>   plt_rbgc%RUPNOB    , &
    RUOH2P     =>   plt_rbgc%RUOH2P    , &
    RUONOB     =>   plt_rbgc%RUONOB    , &
    RUONHB     =>   plt_rbgc%RUONHB    , &
    RUONO3     =>   plt_rbgc%RUONO3    , &
    RUCNHB     =>   plt_rbgc%RUCNHB    , &
    RootNH4Uptake_pvr     =>   plt_rbgc%RootNH4Uptake_pvr    , &
    RootNO3Uptake_pvr     =>   plt_rbgc%RootNO3Uptake_pvr    , &
    RUPNHB     =>   plt_rbgc%RUPNHB    , &
    RCO2N      =>   plt_rbgc%RCO2N     , &
    RootHPO4Uptake_pvr     =>   plt_rbgc%RootHPO4Uptake_pvr    , &
    RUCNH4     =>   plt_rbgc%RUCNH4    , &
    RUOH1B     =>   plt_rbgc%RUOH1B    , &
    RootRespPotential_vr      =>   plt_rbgc%RootRespPotential_vr     , &
    RCO2A      =>   plt_rbgc%RCO2A     , &
    CanPHeight4WatUptake      =>   plt_morph%CanPHeight4WatUptake    , &
    MY         =>   plt_morph%MY       , &
    PrimRootRadius     =>   plt_morph%PrimRootRadius   , &
    PrimRootDepth     =>   plt_morph%PrimRootDepth   , &
    HypoctoylHeight      =>   plt_morph%HypoctoylHeight    , &
    SecndRootRadius      =>   plt_morph%SecndRootRadius    , &
    SecndRootXNum_rpvr      =>   plt_morph%SecndRootXNum_rpvr    , &
    AveSecndRootLen      =>   plt_morph%AveSecndRootLen    , &
    SeedinDepth      =>   plt_morph%SeedinDepth    , &
    NI         =>   plt_morph%NI       , &
    NumRootAxes_pft       =>   plt_morph%NumRootAxes_pft       &
  )

!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER

  RootSinkC_vr=0._R8
  Root1stSink_pvr=0._r8
  Root2ndSink_pvr=0._r8
  RootSinkC=0._r8
  D4995: DO N=1,MY(NZ)
    D4990: DO L=NU,NI(NZ)
!
!     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
!     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
!
!     VLSoilPoreMicP=soil layer volume excluding macropore, rocks
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     RootNH4Uptake_pvr,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RootH2PO4Uptake_pvr,RUPH2B,RootHPO4Uptake_pvr,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!     RUONH4,RUONHB,RUON03,RUONOB=uptake from non-band,band of NH4,NO3 unlimited by O2
!     RUOH2P,RUOH2B,RUOH1P,RUOH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by O2
!     RUCNH4,RUCNHB,RUCN03,RUCNOB=uptake from non-band,band of NH4,NO3 unlimited by nonstructural C
!     RUCH2P,RUCH2B,RUCH1P,RUCH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by nonstructural C
!
      IF(VLSoilPoreMicP(L).GT.ZEROS2)THEN
        CUPRL=0.86_r8*(RootNH4Uptake_pvr(N,L,NZ)+RUPNHB(N,L,NZ) &
          +RootNO3Uptake_pvr(N,L,NZ)+RUPNOB(N,L,NZ)+RootH2PO4Uptake_pvr(N,L,NZ) &
          +RUPH2B(N,L,NZ)+RootHPO4Uptake_pvr(N,L,NZ)+RUPH1B(N,L,NZ))
        CUPRO=0.86_r8*(RUONH4(N,L,NZ)+RUONHB(N,L,NZ) &
          +RUONO3(N,L,NZ)+RUONOB(N,L,NZ)+RUOH2P(N,L,NZ) &
          +RUOH2B(N,L,NZ)+RUOH1P(N,L,NZ)+RUOH1B(N,L,NZ))
        CUPRC=0.86_r8*(RUCNH4(N,L,NZ)+RUCNHB(N,L,NZ) &
          +RUCNO3(N,L,NZ)+RUCNOB(N,L,NZ)+RUCH2P(N,L,NZ) &
          +RUCH2B(N,L,NZ)+RUCH1P(N,L,NZ)+RUCH1B(N,L,NZ))
!
!     ACCUMULATE RESPIRATION IN FLUX ARRAYS
!
!     RCO2A=total root respiration
!     RootRespPotential_vr,RCO2N=RCO2A unltd by O2,nonstructural C
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     CPOOLR=non-structural C mass in root
!
        RootRespPotential_vr(N,L,NZ)=RootRespPotential_vr(N,L,NZ)+CUPRO
        RCO2N(N,L,NZ)=RCO2N(N,L,NZ)+CUPRC
        RCO2A(N,L,NZ)=RCO2A(N,L,NZ)-CUPRL
         RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)-CUPRL
!
!     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
!     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
!     RootNH4Uptake_pvr,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RootH2PO4Uptake_pvr,RUPH2B,RootHPO4Uptake_pvr,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!
        D195: DO K=1,jcplx
          DO NE=1,NumOfPlantChemElmnts
             RootMycoNonstructElmnt_vr(NE,N,L,NZ)= RootMycoNonstructElmnt_vr(NE,N,L,NZ)+RDFOME(NE,N,K,L,NZ)
          ENDDO
        ENDDO D195
         RootMycoNonstructElmnt_vr(ielmn,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmn,N,L,NZ)+RootNH4Uptake_pvr(N,L,NZ)+RUPNHB(N,L,NZ) &
          +RootNO3Uptake_pvr(N,L,NZ)+RUPNOB(N,L,NZ)
         RootMycoNonstructElmnt_vr(ielmp,N,L,NZ)= RootMycoNonstructElmnt_vr(ielmp,N,L,NZ)+RootH2PO4Uptake_pvr(N,L,NZ)+RUPH2B(N,L,NZ) &
          +RootHPO4Uptake_pvr(N,L,NZ)+RUPH1B(N,L,NZ)
!
!     GROWTH OF EACH ROOT AXIS
!
        D4985: DO NR=1,NumRootAxes_pft(NZ)
!
!     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
!
!     RTDP1=primary root depth from soil surface
!     RTDPP=primary root depth from canopy
!     CumSoilThickness=depth from soil surface to layer bottom
!     CanPHeight4WatUptake=canopy height for water uptake
!     RTSK=relative primary root sink strength
!     Root1stSink_pvr=primary root sink strength
!     RootAreaPopu=number of primary root axes
!     RRAD1,SecndRootRadius=primary, secondary root radius
!     RootSinkC,RootSinkC_vr=total root sink strength
!
          IF(N.EQ.ipltroot)THEN
            IF(PrimRootDepth(N,NR,NZ).GT.CumSoilThickness(L-1))THEN
              IF(PrimRootDepth(N,NR,NZ).LE.CumSoilThickness(L))THEN
                RTDPP=PrimRootDepth(N,NR,NZ)+CanPHeight4WatUptake(NZ)
                Root1stSink_pvr(N,L,NR)=RTSK(iPlantMorphologyType(NZ))*RootAreaPopu*PrimRootRadius(N,L,NZ)**2._r8/RTDPP
                RootSinkC(N)=RootSinkC(N)+Root1stSink_pvr(N,L,NR)
                RootSinkC_vr(N,L)=RootSinkC_vr(N,L)+Root1stSink_pvr(N,L,NR)
              ENDIF
            ENDIF
          ENDIF
!
!     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
!     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
!
!     Root1stDepz_vr=depth of primary root axis in layer
!     RTDP1=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     RTDPX=distance behind growing point for secondary roots
!     DLYR=layer thickness
!     SeedinDepth=seeding depth
!     HypoctoylHeight=hypocotyledon height
!     CanPHeight4WatUptake=canopy height for water uptake
!     RTDPS=secondary root depth from canopy
!     RTSKP,RTSKS=primary,secondary root sink strength
!     RTN2=number of secondary root axes
!     Root2ndSink_pvr=total secondary root sink strength
!     AveSecndRootLen=average secondary root length
!     RootSinkC,RootSinkC_vr=total root sink strength
!
          IF(N.EQ.ipltroot)THEN
            Root1stDepz_vr(NR,L)=AZMAX1(PrimRootDepth(ipltroot,NR,NZ)-CumSoilThickness(L-1)-RTDPX)
            Root1stDepz_vr(NR,L)=AZMAX1(AMIN1(DLYR3(L),Root1stDepz_vr(NR,L)) &
              -AZMAX1(SeedinDepth(NZ)-CumSoilThickness(L-1)-HypoctoylHeight(NZ)))
            RTDPS=AMAX1(SeedinDepth(NZ),CumSoilThickness(L-1))+0.5_r8*Root1stDepz_vr(NR,L)+CanPHeight4WatUptake(NZ)
            IF(RTDPS.GT.ZERO)THEN
              RTSKP=RootAreaPopu*PrimRootRadius(N,L,NZ)**2._r8/RTDPS
              RTSKS=safe_adb(SecndRootXNum_rpvr(N,L,NR,NZ)*SecndRootRadius(N,L,NZ)**2._r8,AveSecndRootLen(N,L,NZ))
              IF(RTSKP+RTSKS.GT.ZEROP(NZ))THEN
                Root2ndSink_pvr(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
              ELSE
                Root2ndSink_pvr(N,L,NR)=0._r8
              ENDIF
            ELSE
              Root2ndSink_pvr(N,L,NR)=0._r8
            ENDIF
          ELSE
            Root2ndSink_pvr(N,L,NR)=safe_adb(SecndRootXNum_rpvr(N,L,NR,NZ)*SecndRootRadius(N,L,NZ)**2._r8,AveSecndRootLen(N,L,NZ))
          ENDIF
          RootSinkC(N)=RootSinkC(N)+Root2ndSink_pvr(N,L,NR)
          RootSinkC_vr(N,L)=RootSinkC_vr(N,L)+Root2ndSink_pvr(N,L,NR)
        ENDDO D4985
      ENDIF
    ENDDO D4990
  ENDDO D4995
  end associate
  end subroutine SummarizeRootSink

end module RootMod
