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

  subroutine RootBGCModel(I,J,NZ,BegRemoblize,ICHK1,NRX,PTRT,TFN6_vr,CNRTW,CPRTW,RootAreaPopu)

  implicit none
  integer, intent(in) :: I,J,NZ,BegRemoblize
  integer, intent(inout) :: ICHK1(jroots,JZ1)
  integer, intent(inout)  :: NRX(jroots,JZ1)
  real(r8), intent(in) :: TFN6_vr(JZ1),CNRTW,CPRTW,RootAreaPopu
  real(r8), intent(in):: PTRT      !shoot-root nonstrucal C/N exchange modifier

  integer, parameter  :: NumRootAxes4DeadPlant =0    !
  real(r8) :: TotRootVol
  real(r8) :: fRootGrowPSISense_vr(jroots,JZ1)
  real(r8) :: RootSinkC_vr(jroots,JZ1)
  real(r8) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: RootSinkC(jroots)

  associate(                                &
    SeasonalNonstElms_pft            =>   plt_biom%SeasonalNonstElms_pft      , &
    ZEROL                        =>   plt_biom%ZEROL      , &
    DLYR3                        =>   plt_site%DLYR3      , &
    PlantPopulation_pft          =>   plt_site%PlantPopulation_pft         , &
    ZERO                         =>   plt_site%ZERO       , &
    iPlantPhenolPattern_pft      =>   plt_pheno%iPlantPhenolPattern_pft    , &
    iPlantRootState_pft          =>   plt_pheno%iPlantRootState_pft     , &
    iPlantShootState_pft         =>   plt_pheno%iPlantShootState_pft     , &
    SeedMeanLen_pft              =>   plt_morph%SeedMeanLen_pft      , &
    RootVH2O_pvr                 =>   plt_morph%RootVH2O_pvr    , &
    RootAreaPerPlant_pvr         =>   plt_morph%RootAreaPerPlant_pvr     , &
    RootPoreVol_pvr              =>   plt_morph%RootPoreVol_pvr     , &
    RootLenDensPerPlant_pvr      =>   plt_morph%RootLenDensPerPlant_pvr     , &
    RootLenPerPlant_pvr          =>   plt_morph%RootLenPerPlant_pvr     , &
    NGTopRootLayer_pft           =>   plt_morph%NGTopRootLayer_pft        , &
    RootPorosity_pft             =>   plt_morph%RootPorosity_pft      , &
    SeedVolumeMean_pft           =>   plt_morph%SeedVolumeMean_pft      , &
    SeedAreaMean_pft             =>   plt_morph%SeedAreaMean_pft      , &
    NumRootAxes_pft              =>   plt_morph%NumRootAxes_pft         &
  )
!     ROOT GROWTH
!
  call SummarizeRootSink(NZ,RootAreaPopu,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)

  call RootBiochemistry(I,J,NZ,ICHK1,NRX,TFN6_vr,CNRTW,CPRTW,&
      RootAreaPopu,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,fRootGrowPSISense_vr)
!
!     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
!     GERMINATION)
!
!  if(NZ==1)THEN
!    WRITE(33,*)'ROOTL',RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ),SeedMeanLen_pft(NZ)
!  ELSE
!    WRITE(34,*)'ROOTL',RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ),SeedMeanLen_pft(NZ)    
!  ENDIF

  RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedMeanLen_pft(NZ)
  IF(DLYR3(NGTopRootLayer_pft(NZ)).GT.ZERO)THEN
    !by m3
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)/DLYR3(NGTopRootLayer_pft(NZ))
  ELSE
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=0._r8
  ENDIF
  TotRootVol=RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedVolumeMean_pft(NZ)*PlantPopulation_pft(NZ)
  RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootPorosity_pft(ipltroot,NZ)*TotRootVol
  RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=(1.0_r8-RootPorosity_pft(ipltroot,NZ))*TotRootVol
  RootAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+&
    SeedAreaMean_pft(NZ)

  IF(NumRootAxes4DeadPlant.EQ.NumRootAxes_pft(NZ) .OR. (SeasonalNonstElms_pft(ielmc,NZ).LE.ZEROL(NZ).AND. &
    iPlantPhenolPattern_pft(NZ).NE.iplt_annual))THEN
    iPlantRootState_pft(NZ)=iDead
    iPlantShootState_pft(NZ)=iDead
  ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
  call RootNodulBiochemistry(I,J,NZ,TFN6_vr,fRootGrowPSISense_vr)

  call NonstructlBiomTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,BegRemoblize)
  end associate
  end subroutine RootBGCModel

!------------------------------------------------------------------------------------------

  subroutine RootBiochemistry(I,J,NZ,ICHK1,NRX,TFN6_vr,CNRTW,CPRTW,RootAreaPopu,&
    RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,fRootGrowPSISense_vr)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: ICHK1(jroots,JZ1)
  INTEGER, Intent(inout) :: NRX(jroots,JZ1)
  real(r8), intent(in) :: TFN6_vr(JZ1),CNRTW,CPRTW,RootAreaPopu
  REAL(R8), INTENT(in) :: RootSinkC_vr(jroots,JZ1)
  real(r8), INTENT(in) :: RootSinkC(jroots)
  real(r8), INTENT(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(out) :: fRootGrowPSISense_vr(jroots,JZ1)
  integer :: LL,LZ,L1,L,K,lx,M,NR,N,NTG
  real(r8) :: WFNR
  real(r8) :: WFNRG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: NonstElmGradt
  real(r8) :: CPOOLX
  real(r8) :: DMRTD
  real(r8) :: FWTRT
  real(r8) :: SoilResit4SecndRootPentrate
  real(r8) :: TotRoot2ndLen,TotRoot1stLen
  real(r8) :: TotPopuRoot1stLen_rpvr
  real(r8) :: TotPopuRootLen
  real(r8) :: TotRootVol
  real(r8) :: TotRootArea
  real(r8) :: WTRTTX
  real(r8) :: WTRVCX
  real(r8) :: Root2ndC,Root1stC
  real(r8) :: WTRTLX
  real(r8) :: WTRTTT
  real(r8) :: WTRTT
  real(r8) :: XFRC

!     begin_execution
  associate(                            &
    RootMycoNonstElms_rpvr          =>   plt_biom%RootMycoNonstElms_rpvr     , &
    RootMycoActiveBiomC_pvr       =>   plt_biom%RootMycoActiveBiomC_pvr     , &
    RootElms_pft                  =>   plt_biom%RootElms_pft      , &
    ZEROP                         =>   plt_biom%ZEROP      , &
    SeasonalNonstElms_pft             =>   plt_biom%SeasonalNonstElms_pft      , &
    FWODRE                        =>   plt_allom%FWODRE    , &
    RootBiomGrowthYield           =>   plt_allom%RootBiomGrowthYield     , &
    iPlantRootProfile_pft         =>   plt_pheno%iPlantRootProfile_pft    , &
    PSIRoot_pvr                   =>   plt_ew%PSIRoot_pvr       , &
    PSIRootTurg_vr                =>   plt_ew%PSIRootTurg_vr        , &
    RootGasLossDisturb_pft        =>   plt_bgcr%RootGasLossDisturb_pft    , &
    trcg_rootml_pvr               =>   plt_rbgc%trcg_rootml_pvr       , &
    trcs_rootml_pvr               =>   plt_rbgc%trcs_rootml_pvr, &
    SoilResit4RootPentrate_vr     =>   plt_soilchem%SoilResit4RootPentrate_vr   , &
    VLSoilPoreMicP                =>   plt_soilchem%VLSoilPoreMicP   , &
    NU                            =>   plt_site%NU         , &
    ZERO                          =>   plt_site%ZERO       , &
    PlantPopulation_pft           =>   plt_site%PlantPopulation_pft         , &
    ZEROS2                        =>   plt_site%ZEROS2     , &
    DLYR3                         =>   plt_site%DLYR3      , &
    NL                            =>   plt_site%NL         , &
    k_fine_litr                   =>   pltpar%k_fine_litr  , &
    RootPoreVol_pvr               =>   plt_morph%RootPoreVol_pvr     , &
    RootLenDensPerPlant_pvr       =>   plt_morph%RootLenDensPerPlant_pvr     , &
    RootAreaPerPlant_pvr          =>   plt_morph%RootAreaPerPlant_pvr     , &
    RootVH2O_pvr                  =>   plt_morph%RootVH2O_pvr    , &
    Max1stRootRadius_pft1         =>   plt_morph%Max1stRootRadius_pft1    , &
    Max2ndRootRadius_pft1         =>   plt_morph%Max2ndRootRadius_pft1    , &
    Max1stRootRadius_pft          =>   plt_morph%Max1stRootRadius_pft    , &
    Root1stRadius_pvr             =>   plt_morph%Root1stRadius_pvr    , &
    RootLenPerPlant_pvr           =>   plt_morph%RootLenPerPlant_pvr     , &
    AveLen2ndRoot_pvr             =>   plt_morph%AveLen2ndRoot_pvr     , &
    Radius2ndRoot_pvr             =>   plt_morph%Radius2ndRoot_pvr     , &
    Root2ndXSecArea_pft           =>   plt_morph%Root2ndXSecArea_pft    , &
    Max2ndRootRadius_pft          =>   plt_morph%Max2ndRootRadius_pft    , &
    Root2ndXNum_pvr               =>   plt_morph%Root2ndXNum_pvr      , &
    Root1stXSecArea_pft           =>   plt_morph%Root1stXSecArea_pft    , &
    RootPorosity_pft              =>   plt_morph%RootPorosity_pft     , &
    MaxSoiL4Root                  =>   plt_morph%MaxSoiL4Root        , &
    RootVolPerMassC_pft           =>   plt_morph%RootVolPerMassC_pft      , &
    SeedMeanLen_pft               =>   plt_morph%SeedMeanLen_pft      , &
    MY                            =>   plt_morph%MY        , &
    NGTopRootLayer_pft            =>   plt_morph%NGTopRootLayer_pft        , &
    NIXBotRootLayer_pft           =>   plt_morph%NIXBotRootLayer_pft        &
  )

  NIXBotRootLayer_pft(NZ)=NGTopRootLayer_pft(NZ)
!
!
!     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
!
  D5010: DO N=1,MY(NZ)
    D5000: DO L=NU,MaxSoiL4Root(NZ)
!
!     IDENTIFY NEXT LOWER ROOT LAYER
!
!     VLSoilPoreMicP=soil layer volume excluding macropore, rocks
!
      IF(VLSoilPoreMicP(L).GT.ZEROS2)THEN
        !why not using MaxSoiL4Root(NZ)
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
!     SoilResit4RootPentrate_vr,SoilResit4SecndRootPentrate=soil resistance to secondary root penetration (MPa)
!     Radius2ndRoot_pvr=secondary root radius
!     WFNR=water function for root extension
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     fRootGrowPSISense_vr,WFNRG=growth,respiration function of root water potential
!     PSIRoot_pvr,PSIRootTurg_vr=root total,turgor water potential
!     DMRT=root growth yield
!
        SoilResit4SecndRootPentrate=SoilResit4RootPentrate_vr(L)*Radius2ndRoot_pvr(N,L,NZ)/1.0E-03_r8
        WFNR=AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-PSIMin4OrganExtens-SoilResit4SecndRootPentrate))
        IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
          fRootGrowPSISense_vr(N,L)=EXP(0.05_r8*PSIRoot_pvr(N,L,NZ))
          WFNRG=WFNR**0.10_r8
        ELSE
          fRootGrowPSISense_vr(N,L)=EXP(0.10_r8*PSIRoot_pvr(N,L,NZ))
          WFNRG=WFNR**0.25_r8
        ENDIF
        DMRTD=1.0_r8-RootBiomGrowthYield(NZ)
!
!     FOR EACH ROOT AXIS
!
        call GrowRootMycoAxes(N,L,L1,NZ,NRX,fRootGrowPSISense_vr,ICHK1,WFNR,WFNRG,TFN6_vr,RootAreaPopu,DMRTD,&
          RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,CNRTW,CPRTW,TotRoot1stLen,Root2ndC,Root1stC,&
          TotRoot2ndLen)

!
!     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
!     SEASONAL STORAGE POOL IS DEPLETED
!
!     RootMycoActiveBiomC_pvr,WTRT=total root C mass
!     WTRVC=storage C
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!     CPOOLR=non-structural C mass in root
!
        IF(L.LE.NIXBotRootLayer_pft(NZ))THEN
          IF(RootMycoActiveBiomC_pvr(N,L,NZ).GT.ZEROP(NZ).AND. &
            RootElms_pft(ielmc,NZ).GT.ZEROP(NZ) &
            .AND.SeasonalNonstElms_pft(ielmc,NZ).LT.XFRX*RootElms_pft(ielmc,NZ))THEN
            FWTRT=RootMycoActiveBiomC_pvr(N,L,NZ)/RootElms_pft(ielmc,NZ)
            WTRTLX=RootMycoActiveBiomC_pvr(N,L,NZ)
            WTRTTX=RootElms_pft(ielmc,NZ)*FWTRT
            WTRTTT=WTRTLX+WTRTTX
            CPOOLX=AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
            WTRVCX=AZMAX1(SeasonalNonstElms_pft(ielmc,NZ)*FWTRT)
            NonstElmGradt=(WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
            XFRC=AZMIN1(XFRY*NonstElmGradt)
            RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+XFRC
            SeasonalNonstElms_pft(ielmc,NZ)=SeasonalNonstElms_pft(ielmc,NZ)-XFRC
          ENDIF
        ENDIF
!
!     ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
!     TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
!
!     TotRoot1stLen=total primary root length
!     Root1stC=total primary root C mass
!     TotRoot2ndLen=total secondary root length
!     Root2ndC=total secondary root C mass
!     TotPopuRootLen=total root length
!     WTRTT=total root C mass
!     FWOOD=C woody fraction in root:0=woody,1=non-woody
!     PP=PFT population
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     TotRootVol,RootVH2O_pvr,RootPoreVol_pvr=root or myco total,aqueous,gaseous volume
!     RRAD1,Radius2ndRoot_pvr=primary,secondary root radius
!     RootAreaPerPlant_pvr=root surface area per plant
!     AveLen2ndRoot_pvr=average secondary root length
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!
        IF(N.EQ.ipltroot)THEN
          TotRoot1stLen=TotRoot1stLen*FWODRE(ielmc,k_fine_litr)
          TotRoot2ndLen=TotRoot2ndLen*FWODRE(ielmc,k_fine_litr)
        ENDIF
        TotPopuRoot1stLen_rpvr=TotRoot1stLen*PlantPopulation_pft(NZ)
        TotPopuRootLen=TotRoot2ndLen+TotPopuRoot1stLen_rpvr
        WTRTT=Root2ndC+Root1stC
        IF(TotPopuRootLen.GT.ZEROP(NZ) .AND. WTRTT.GT.ZEROP(NZ) .AND. PlantPopulation_pft(NZ).GT.ZEROP(NZ))THEN
          RootLenPerPlant_pvr(N,L,NZ)=TotPopuRootLen/PlantPopulation_pft(NZ)
          IF(DLYR3(L).GT.ZERO)THEN
            !per volume
            RootLenDensPerPlant_pvr(N,L,NZ)=RootLenPerPlant_pvr(N,L,NZ)/DLYR3(L)
          ELSE
            RootLenDensPerPlant_pvr(N,L,NZ)=0._r8
          ENDIF
          TotRootVol=AMAX1(Root1stXSecArea_pft(N,NZ)*TotPopuRoot1stLen_rpvr+Root2ndXSecArea_pft(N,NZ)*TotRoot2ndLen &
            ,WTRTT*RootVolPerMassC_pft(N,NZ)*PSIRootTurg_vr(N,L,NZ))
          RootPoreVol_pvr(N,L,NZ)=RootPorosity_pft(N,NZ)*TotRootVol
          RootVH2O_pvr(N,L,NZ)=(1.0_r8-RootPorosity_pft(N,NZ))*TotRootVol
          !primary roots
          Root1stRadius_pvr(N,L,NZ)=AMAX1(Max1stRootRadius_pft1(N,NZ),&
            (1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Max1stRootRadius_pft(N,NZ))
          !secondary roots
          Radius2ndRoot_pvr(N,L,NZ)=AMAX1(Max2ndRootRadius_pft1(N,NZ),&
            (1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Max2ndRootRadius_pft(N,NZ))
          TotRootArea=TwoPiCON*(Root1stRadius_pvr(N,L,NZ)*TotPopuRoot1stLen_rpvr+Radius2ndRoot_pvr(N,L,NZ)*TotRoot2ndLen)
          IF(Root2ndXNum_pvr(N,L,NZ).GT.ZEROP(NZ))THEN
            AveLen2ndRoot_pvr(N,L,NZ)=AMAX1(MinAve2ndRootLen,TotRoot2ndLen/Root2ndXNum_pvr(N,L,NZ))
          ELSE
            AveLen2ndRoot_pvr(N,L,NZ)=MinAve2ndRootLen
          ENDIF
          RootAreaPerPlant_pvr(N,L,NZ)=TotRootArea/PlantPopulation_pft(NZ)

        ELSE
          RootLenPerPlant_pvr(N,L,NZ)=0._r8
          RootLenDensPerPlant_pvr(N,L,NZ)=0._r8
          RootPoreVol_pvr(N,L,NZ)=0._r8
          RootVH2O_pvr(N,L,NZ)=0._r8
          Root1stRadius_pvr(N,L,NZ)=Max1stRootRadius_pft(N,NZ)
          Radius2ndRoot_pvr(N,L,NZ)=Max2ndRootRadius_pft(N,NZ)
          RootAreaPerPlant_pvr(N,L,NZ)=0._r8
          AveLen2ndRoot_pvr(N,L,NZ)=MinAve2ndRootLen
          DO NTG=idg_beg,idg_end-1
            RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ) &
              -(trcg_rootml_pvr(NTG,N,L,NZ)+trcs_rootml_pvr(NTG,N,L,NZ))
          ENDDO
          trcg_rootml_pvr(idg_beg:idg_end-1,N,L,NZ)=0._r8
          trcs_rootml_pvr(idg_beg:idg_end-1,N,L,NZ)=0._r8
        ENDIF
      ENDIF
    ENDDO D5000
  ENDDO D5010
  end associate
  end subroutine RootBiochemistry

!------------------------------------------------------------------------------------------

  subroutine GrowRootMycoAxes(N,L,L1,NZ,NRX,fRootGrowPSISense_vr,ICHK1,WFNR,WFNRG,TFN6_vr,&
    RootAreaPopu,DMRTD,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,CNRTW,CPRTW,&
    TotRoot2ndLen,TotRoot1stLen,Root2ndC,Root1stC)
  implicit none
  INTEGER, INTENT(IN) :: N,L,L1,NZ
  integer, intent(inout) :: NRX(jroots,JZ1)
  real(r8), intent(in) :: TFN6_vr(JZ1),RootAreaPopu
  real(r8), intent(in) :: fRootGrowPSISense_vr(jroots,JZ1)
  real(r8), intent(in) :: DMRTD
  REAL(R8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: CNRTW,CPRTW
  integer, intent(inout) :: ICHK1(jroots,JZ1)
  real(r8), intent(inout):: WFNR,WFNRG
  real(r8), intent(out) :: TotRoot2ndLen,TotRoot1stLen,Root2ndC,Root1stC
  real(r8) :: CNRDA,CNRDM
  real(r8) :: CNPG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: FSNC2
  real(r8) :: ZADD1M,RootMycoNonst4Grow(NumPlantChemElms),ZADD2M
  real(r8) :: ZPOOLB
  real(r8) :: RootMycoNonst4Grow_Oulm
  real(r8) :: DMRTR
  real(r8) :: FNP
  real(r8) :: FRTN
  real(r8) :: FRCO2
  real(r8) :: FSNCM
  real(r8) :: FSNCP
  real(r8) :: GRTWGM
  real(r8) :: Root2ndExtension
  real(r8) :: RootCYldO2ltd
  real(r8) :: RootNetGrowthElms(NumPlantChemElms)
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y
  real(r8) :: RGrowCO2_Oltd,RCO2T_Oltd
  real(r8) :: RCO2XM
  real(r8) :: RCO2YM
  real(r8) :: RGrowCO2_OUlm
  real(r8) :: RCO2T_Oulm
  real(r8) :: Rmaint_CO2,RNonstCO2_OUlm,RNonstCO2_Oltd
  real(r8) :: RootChemElmRemob(NumPlantChemElms)
  real(r8) :: RTN2X,RTN2Y
  REAL(R8) :: SNCR,SNCRM
  real(r8) :: TFRCO2
  real(r8) :: RCCC,RCCN,RCCP
  integer :: NR,M,NE

!begin_execution
  associate(                          &
    RootMyco2ndStrutElms_rpvr            =>  plt_biom%RootMyco2ndStrutElms_rpvr     , &
    RootNonstructElmConc_pvr        =>  plt_biom%RootNonstructElmConc_pvr     , &
    RootMyco1stStrutElms_rpvr           =>  plt_biom%RootMyco1stStrutElms_rpvr     , &
    RootMycoNonstElms_rpvr            =>  plt_biom%RootMycoNonstElms_rpvr     , &
    RootProteinC_pvr                =>  plt_biom%RootProteinC_pvr      , &
    RootMycoActiveBiomC_pvr         =>  plt_biom%RootMycoActiveBiomC_pvr     , &
    ZEROP                           =>  plt_biom%ZEROP      , &
    RCO2A_pvr                       =>  plt_rbgc%RCO2A_pvr      , &
    RCO2N_pvr                       =>  plt_rbgc%RCO2N_pvr      , &
    RootRespPotent_pvr              =>  plt_rbgc%RootRespPotent_pvr      , &
    RAutoRootO2Limter_pvr          =>  plt_rbgc%RAutoRootO2Limter_pvr       , &
    LitfalStrutElms_pvr               =>  plt_bgcr%LitfalStrutElms_pvr       , &
    rCNNonstructRemob_pft           =>  plt_allom%rCNNonstructRemob_pft     , &
    rCPNonstructRemob_pft           =>  plt_allom%rCPNonstructRemob_pft      , &
    FWODRE                          =>  plt_allom%FWODRE    , &
    CNRTS                           =>  plt_allom%CNRTS     , &
    CPRTS                           =>  plt_allom%CPRTS     , &
    RootBiomGrowthYield             =>  plt_allom%RootBiomGrowthYield      , &
    k_woody_litr                    =>  pltpar%k_woody_litr,&
    k_fine_litr                     =>  pltpar%k_fine_litr, &
    icwood                          =>  pltpar%icwood       , &
    iroot                           =>  pltpar%iroot        , &
    inonstruct                      =>  pltpar%inonstruct     , &
    iPlantRootProfile_pft           =>  plt_pheno%iPlantRootProfile_pft    , &
    iPlantPhenolType_pft            =>  plt_pheno%iPlantPhenolType_pft    , &
    fTgrowRootP_vr                  =>  plt_pheno%fTgrowRootP_vr      , &
    iPlantCalendar_brch             =>  plt_pheno%iPlantCalendar_brch    , &
    CFOPE                           =>  plt_soilchem%CFOPE  , &
    SoilResit4RootPentrate_vr       =>  plt_soilchem%SoilResit4RootPentrate_vr   , &
    DLYR3                           =>  plt_site%DLYR3      , &
    ZERO                            =>  plt_site%ZERO       , &
    Root2ndXNum_pvr                 =>  plt_morph%Root2ndXNum_pvr      , &
    MaxSeedCMass                    =>  plt_morph%MaxSeedCMass      , &
    NGTopRootLayer_pft              =>  plt_morph%NGTopRootLayer_pft       , &
    NIXBotRootLayer_pft             =>  plt_morph%NIXBotRootLayer_pft      , &
    Root2ndXNum_rpvr                =>  plt_morph%Root2ndXNum_rpvr     , &
    NumRootAxes_pft                 =>  plt_morph%NumRootAxes_pft      , &
    RootBranchFreq_pft              =>  plt_morph%RootBranchFreq_pft      , &
    SeedDepth_pft                   =>  plt_morph%SeedDepth_pft     , &
    Root1stLen_rpvr                 =>  plt_morph%Root1stLen_rpvr     , &
    Root2ndSpecLen_pft              =>  plt_morph%Root2ndSpecLen_pft    , &
    Root2ndLen_pvr                  =>  plt_morph%Root2ndLen_pvr     , &
    NIXBotRootLayer_rpft            =>  plt_morph%NIXBotRootLayer_rpft      , &
    MainBranchNum_pft               =>  plt_morph%MainBranchNum_pft       , &
    C4PhotosynDowreg_brch           =>  plt_photo%C4PhotosynDowreg_brch       &
  )
  TotRoot2ndLen=0._r8
  TotRoot1stLen=0._r8
  Root2ndC=0._r8
  Root1stC=0._r8
  D5050: DO NR=1,NumRootAxes_pft(NZ)
!
!     SECONDARY ROOT EXTENSION
!   make sure current layer is no deeper than the lowest root layer
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
      IF(RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
        CNPG=AMIN1(RootNonstructElmConc_pvr(ielmn,N,L,NZ)/(RootNonstructElmConc_pvr(ielmn,N,L,NZ) &
          +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI),RootNonstructElmConc_pvr(ielmp,N,L,NZ) &
          /(RootNonstructElmConc_pvr(ielmp,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI))
      ELSE
        CNPG=1.0_r8
      ENDIF
!
!     SECONDARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     Rmaint_CO2=root maintenance respiration
!     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT2N=secondary root N mass
!     TFN6_vr=temperature function for root maintenance respiration
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     fRootGrowPSISense_vr=growth function of root water potential
!
      Rmaint_CO2=AZMAX1(RmSpecPlant*RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ))*TFN6_vr(L)
      IF(is_root_shallow(iPlantRootProfile_pft(NZ)).OR. &
        iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
        Rmaint_CO2=Rmaint_CO2*fRootGrowPSISense_vr(N,L)
      ENDIF
!
!     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RNonstCO2_OUlm=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     
!     FRTN=fraction of secondary root sink strength in axis
!     CPOOL=non-structural C mass
!     fTgrowRootP_vr=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!     fRootGrowPSISense_vr=growth function of root water potential
!
      RNonstCO2_OUlm=AZMAX1(VMXC*FRTN*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)) &
        *fTgrowRootP_vr(L,NZ)*fRootGrowPSISense_vr(N,L)&
        *CNPG*C4PhotosynDowreg_brch(MainBranchNum_pft(NZ),NZ)
!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RAutoRootO2Limter_pvr=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RNonstCO2_Oltd=RNonstCO2_OUlm*RAutoRootO2Limter_pvr(N,L,NZ)
      RCO2XM=RNonstCO2_OUlm-Rmaint_CO2
      RCO2X=RNonstCO2_Oltd-Rmaint_CO2
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
!     RGrowCO2_OUlm,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))
      PPOOLB=AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ),PPOOLB*DMRTR/CPRTS(NZ))
      IF(RCO2YM.GT.0.0_r8)THEN
        RGrowCO2_OUlm=AMIN1(RCO2YM,FNP)
      ELSE
        RGrowCO2_OUlm=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
        RGrowCO2_Oltd=AMIN1(RCO2Y,FNP*RAutoRootO2Limter_pvr(N,L,NZ))
      ELSE
        RGrowCO2_Oltd=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN 'READQ'
!
!     RootMycoNonst4Grow_Oulm,RootMycoNonst4Grow(ielmc)=total non-structural C used in growth and respn unltd,ltd by O2
!     RGrowCO2_OUlm,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,RootCYldO2ltd=root C growth unltd,ltd by O2
!     RootBiomGrowthYield=root growth yield
!     ZADD2M,RootMycoNonst4Grow(ielmn),RootMycoNonst4Grow(ielmp)=nonstructural N,P unlimited,limited by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
      RootMycoNonst4Grow_Oulm=RGrowCO2_OUlm/DMRTD
      RootMycoNonst4Grow(ielmc)=RGrowCO2_Oltd/DMRTD
      GRTWGM=RootMycoNonst4Grow_Oulm*RootBiomGrowthYield(NZ)
      RootCYldO2ltd=RootMycoNonst4Grow(ielmc)*RootBiomGrowthYield(NZ)
      ZADD2M=AZMAX1(GRTWGM*CNRTW)
      RootMycoNonst4Grow(ielmn)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootCYldO2ltd*CNRTW))
      RootMycoNonst4Grow(ielmp)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootCYldO2ltd*CPRTW))
      CNRDM=AZMAX1(1.70_r8*ZADD2M)
      CNRDA=AZMAX1(1.70_r8*RootMycoNonst4Grow(ielmn))
!
!     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LitrFall
!
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
      IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0.AND.&
        RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
        CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmn,N,L,NZ)&
          ,RootNonstructElmConc_pvr(ielmn,N,L,NZ) &
          +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI) &
          ,safe_adb(RootNonstructElmConc_pvr(ielmp,N,L,NZ),RootNonstructElmConc_pvr(ielmp,N,L,NZ) &
          +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI)))
        CNC=AZMAX1(AMIN1(1.0_r8 &
          ,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ),RootNonstructElmConc_pvr(ielmc,N,L,NZ)+ &
          RootNonstructElmConc_pvr(ielmn,N,L,NZ)/CNKI)))
        CPC=AZMAX1(AMIN1(1.0_r8 &
          ,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ),RootNonstructElmConc_pvr(ielmc,N,L,NZ)+ &
          RootNonstructElmConc_pvr(ielmp,N,L,NZ)/CPKI)))
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
!     RAutoRootO2Limter_pvr=constraint by O2 consumption on all root processes
!     RootChemElmRemob(ielmc),RootChemElmRemob(ielmn),RootChemElmRemob(ielmp)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC2=fraction of secondary root C to be remobilized
!
      IF(-RCO2XM.GT.0.0_r8)THEN
        IF(-RCO2XM.LT.RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
          SNCRM=-RCO2XM
        ELSE
          SNCRM=AZMAX1(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)
        ENDIF
      ELSE
        SNCRM=0._r8
      ENDIF
      IF(-RCO2X.GT.0.0_r8)THEN
        IF(-RCO2X.LT.RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
          SNCR=-RCO2X
        ELSE
          SNCR=AZMAX1(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)*RAutoRootO2Limter_pvr(N,L,NZ)
        ENDIF
      ELSE
        SNCR=0._r8
      ENDIF

      IF(SNCR.GT.0.0_r8 .AND. RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZEROP(NZ))THEN
        RootChemElmRemob(ielmc)=RCCC*RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
        RootChemElmRemob(ielmn)=RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ)*(RCCN+(1.0_r8-RCCN)* &
          RootChemElmRemob(ielmc)/RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ))
        RootChemElmRemob(ielmp)=RootMyco2ndStrutElms_rpvr(ielmp,N,L,NR,NZ)*(RCCP+(1.0_r8-RCCP)* &
          RootChemElmRemob(ielmc)/RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ))
        IF(RootChemElmRemob(ielmc).GT.ZEROP(NZ))THEN
          FSNC2=AZMAX1(AMIN1(1.0_r8,SNCR/RootChemElmRemob(ielmc)))
        ELSE
          FSNC2=1.0_r8
        ENDIF
      ELSE
        RootChemElmRemob(1:NumPlantChemElms)=0._r8
        FSNC2=0._r8
      ENDIF
!
!     SECONDARY ROOT LitrFall CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC2=fraction of secondary root C to be remobilized
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RootChemElmRemob(ielmc),RootChemElmRemob(ielmn),RootChemElmRemob(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!      
      D6350: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ) &
            +CFOPE(NE,icwood,M,NZ) &
            *FSNC2*(RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)-RootChemElmRemob(NE))*FWODRE(NE,k_woody_litr)

          LitfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
            +CFOPE(NE,iroot,M,NZ) &
            *FSNC2*(RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)-RootChemElmRemob(NE))*FWODRE(NE,k_fine_litr)
        ENDDO    
      ENDDO D6350
      
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     Rmaint_CO2=root maintenance respiration
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RootMycoNonst4Grow(ielmc)=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC2=fraction of secondary root C to be remobilized
!     RootChemElmRemob(ielmc),RootChemElmRemob(ielmn),RootChemElmRemob(ielmp)=remobilization of C,N,P from senescing root
!     RootMycoNonst4Grow(ielmn),RootMycoNonst4Grow(ielmp)=nonstructural N,P ltd by O2 used in growth
!
       RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-AMIN1(Rmaint_CO2,RNonstCO2_Oltd) &
        -CNRDA-SNCR

       DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ) &
          -RootMycoNonst4Grow(NE)+FSNC2*RootChemElmRemob(NE)
       ENDDO
!
!     TOTAL SECONDARY ROOT RESPIRATION
!
!     RCO2T_Oulm,RCO2T_Oltd=total C respiration unltd,ltd by O2
!     Rmaint_CO2=root maintenance respiration
!     RNonstCO2_OUlm,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
!     RGrowCO2_OUlm,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A_pvr=total root respiration
!     RootRespPotent_pvr,RCO2N_pvr=RCO2A_pvr unltd by O2,nonstructural C
!
      RCO2T_Oulm=AMIN1(Rmaint_CO2,RNonstCO2_OUlm)+RGrowCO2_OUlm+SNCRM+CNRDM
      RCO2T_Oltd=AMIN1(Rmaint_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+SNCR+CNRDA
      RootRespPotent_pvr(N,L,NZ)=RootRespPotent_pvr(N,L,NZ)+RCO2T_Oulm
      RCO2N_pvr(N,L,NZ)=RCO2N_pvr(N,L,NZ)+RCO2T_Oltd
      RCO2A_pvr(N,L,NZ)=RCO2A_pvr(N,L,NZ)-RCO2T_Oltd
!
!     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     Root2ndExtension=secondary root length extension
!     RootCYldO2ltd=secondary root C growth ltd by O2
!     Root2ndSpecLen_pft=specific secondary root length from startq.f
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     FSNC2=fraction of secondary root C to be remobilized
!     Root2ndLen_pvr=secondary root length
!     RootNetGrowthElms(ielmc),RootNetGrowthElms(ielmn),RootNetGrowthElms(ielmp)=net root C,N,P growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RootMycoNonst4Grow(ielmn),RootMycoNonst4Grow(ielmp)=nonstructural N,P ltd by O2 used in growth
!
      Root2ndExtension=RootCYldO2ltd*Root2ndSpecLen_pft(N,NZ)*WFNR*FWODRE(ielmc,k_fine_litr) &
        -FSNC2*Root2ndLen_pvr(N,L,NR,NZ)
      RootNetGrowthElms(ielmc)=RootCYldO2ltd-FSNC2*RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
      RootNetGrowthElms(ielmn)=RootMycoNonst4Grow(ielmn)-FSNC2*RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ)
      RootNetGrowthElms(ielmp)=RootMycoNonst4Grow(ielmp)-FSNC2*RootMyco2ndStrutElms_rpvr(ielmp,N,L,NR,NZ)
!
!     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
!     AND AXIS NUMBER
!
!     Root2ndLen_pvr=secondary root length
!     Root2ndExtension=secondary root length extension
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RootNetGrowthElms(ielmc),RootNetGrowthElms(ielmn),RootNetGrowthElms(ielmp)=net root C,N,P growth
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!     RootBranchFreq_pft=root branching frequency from PFT file
!     Root2ndXNum_rpvr,Root2ndXNum_pvr=number of secondary root axes
!
      Root2ndLen_pvr(N,L,NR,NZ)=Root2ndLen_pvr(N,L,NR,NZ)+Root2ndExtension
      DO NE=1,NumPlantChemElms
        RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)&
          +RootNetGrowthElms(NE)
      ENDDO
      RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)+AMIN1(rCNNonstructRemob_pft(NZ)&
        *RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ) &
        ,rCPNonstructRemob_pft(NZ)*RootMyco2ndStrutElms_rpvr(ielmp,N,L,NR,NZ))
      TotRoot2ndLen=TotRoot2ndLen+Root2ndLen_pvr(N,L,NR,NZ)
      Root2ndC=Root2ndC+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)

      RTN2X=RootBranchFreq_pft(NZ)*RootAreaPopu
      RTN2Y=RootBranchFreq_pft(NZ)*RTN2X
      Root2ndXNum_rpvr(N,L,NR,NZ)=(RTN2X+RTN2Y)*DLYR3(L)
      Root2ndXNum_pvr(N,L,NZ)=Root2ndXNum_pvr(N,L,NZ)+Root2ndXNum_rpvr(N,L,NR,NZ)
!
!     PRIMARY ROOT EXTENSION
!
!     SoiBulkDensity=soil bulk density
!     RTDP1,Root1stDepz2Surf=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     ICHKL=flag for identifying layer with primary root tip
!     RTN1=number of primary root axes
!     RootAreaPopu=multiplier for number of primary root axes
!
      IF(N.EQ.ipltroot)THEN
        call GrowRootAxes(N,NR,L,L1,NZ,RootAreaPopu,Root1stSink_pvr,&
          Root2ndSink_pvr,RootSinkC_vr,RootChemElmRemob,fRootGrowPSISense_vr,TFN6_vr,&
          SNCRM,SNCR,DMRTD,CNRTW,CPRTW,TotRoot1stLen,Root1stC,ICHK1,NRX)
      ENDIF
      TotRoot1stLen=TotRoot1stLen+Root1stLen_rpvr(N,L,NR,NZ)
      Root1stC=Root1stC+RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
!     ENDIF
    ENDIF
    NIXBotRootLayer_pft(NZ)=MAX(NIXBotRootLayer_pft(NZ),NIXBotRootLayer_rpft(NR,NZ))
  ENDDO D5050
  end associate
  end subroutine GrowRootMycoAxes

!------------------------------------------------------------------------------------------

  subroutine GrowRootAxes(N,NR,L,L1,NZ,RootAreaPopu,Root1stSink_pvr,&
    Root2ndSink_pvr,RootSinkC_vr,RootChemElmRemob,fRootGrowPSISense_vr,TFN6_vr,&
    SNCRM,SNCR,DMRTD,CNRTW,CPRTW,TotRoot1stLen,Root1stC,ICHK1,NRX)
  !
  !grow root axes
  implicit none
  integer,  intent(in) :: N,NR,L,L1,NZ
  real(r8), intent(in) :: RootAreaPopu
  REAL(R8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense_vr(jroots,JZ1)    
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: RootChemElmRemob(NumPlantChemElms)  
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: SNCRM,SNCR  
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRTD  
  integer,  intent(inout) :: ICHK1(jroots,JZ1)
  integer,  intent(inout) :: NRX(jroots,JZ1)  
  real(r8), intent(inout) :: TotRoot1stLen
  real(r8), intent(inout) :: Root1stC
  real(r8) :: Root1stDepz2Surf
  real(r8) :: FRTN
  real(r8) :: FSNC1  
  real(r8) :: GRTWTM  
  real(r8) :: WFNR,WFNRG
  real(r8) :: SoilResit4PrimRootPentration  
  real(r8) :: CNPG,TFRCO2,FRCO2  
  real(r8) :: Rmaint_CO2,RootCYldO2ltd
  real(r8) :: RNonstCO2_OUlm,RGrowCO2_Oltd,RGrowCO2_OUlm
  real(r8) :: RNonstCO2_Oltd,RCO2XM,RCO2X,RCO2YM,RCO2Y
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4Grow_Oulm,RootMycoNonst4Grow(NumPlantChemElms)
  real(r8) :: GRTWGM,ZADD1M,CNRDM,CNRDA
  real(r8) :: RCO2T_Oulm,RCO2T_Oltd
  real(r8) :: FSNCM,FSNCP  
  real(r8) :: RootNetGrowthElms(NumPlantChemElms)
  integer  :: NE,LL,LX,M
  
  associate(                                                            &
    ZERO                            =>  plt_site%ZERO       , &  
    RCO2A_pvr                       =>  plt_rbgc%RCO2A_pvr      , &
    RCO2N_pvr                       =>  plt_rbgc%RCO2N_pvr      , &    
    RootRespPotent_pvr              =>  plt_rbgc%RootRespPotent_pvr      , &        
    LitfalStrutElms_pvr               =>  plt_bgcr%LitfalStrutElms_pvr       , &    
    iPlantRootProfile_pft           =>  plt_pheno%iPlantRootProfile_pft    , &   
    iPlantPhenolType_pft            =>  plt_pheno%iPlantPhenolType_pft    , & 
    fTgrowRootP_vr                  =>  plt_pheno%fTgrowRootP_vr      , &    
    RAutoRootO2Limter_pvr           =>  plt_rbgc%RAutoRootO2Limter_pvr       , &     
    MaxNumRootLays                  =>  plt_site%MaxNumRootLays         , &
    CumSoilThickness                =>  plt_site%CumSoilThickness     , &         
    RootMycoNonstElms_rpvr            =>  plt_biom%RootMycoNonstElms_rpvr     , &    
    Root1stElm_raxs                 =>  plt_biom%Root1stElm_raxs     , &    
    ZEROP                           =>  plt_biom%ZEROP      , &        
    RootMycoActiveBiomC_pvr         =>  plt_biom%RootMycoActiveBiomC_pvr     , &    
    RootNonstructElmConc_pvr        =>  plt_biom%RootNonstructElmConc_pvr     , &        
    RootMyco2ndStrutElms_rpvr            =>  plt_biom%RootMyco2ndStrutElms_rpvr     , &  
    RootMyco1stStrutElms_rpvr           =>  plt_biom%RootMyco1stStrutElms_rpvr      , &    
    icwood                          =>  pltpar%icwood       , &
    iroot                           =>  pltpar%iroot        , &    
    inonstruct                      =>  pltpar%inonstruct     , &    
    k_woody_litr                    =>  pltpar%k_woody_litr,&
    k_fine_litr                     =>  pltpar%k_fine_litr, &      
    CNRTS                           =>  plt_allom%CNRTS     , &
    CPRTS                           =>  plt_allom%CPRTS     , &    
    FWODRE                          =>  plt_allom%FWODRE    , &    
    RootBiomGrowthYield             =>  plt_allom%RootBiomGrowthYield      , &    
    Root1stLen_rpvr                 =>  plt_morph%Root1stLen_rpvr     , &    
    SeedDepth_pft                   =>  plt_morph%SeedDepth_pft     , &    
    NGTopRootLayer_pft              =>  plt_morph%NGTopRootLayer_pft       , &    
    Root2ndLen_pvr                  =>  plt_morph%Root2ndLen_pvr     , &    
    Root1stXNumL_pvr                =>  plt_morph%Root1stXNumL_pvr     , &    
    Root1stRadius_pvr               =>  plt_morph%Root1stRadius_pvr    , &    
    Root1stDepz_pft                 =>  plt_morph%Root1stDepz_pft    , &    
    MainBranchNum_pft               =>  plt_morph%MainBranchNum_pft       , &
    NIXBotRootLayer_rpft            =>  plt_morph%NIXBotRootLayer_rpft      , &     
    PSIRootTurg_vr                  =>  plt_ew%PSIRootTurg_vr        , &       
    CFOPE                           =>  plt_soilchem%CFOPE  , &    
    SoiBulkDensity                  =>  plt_soilchem%SoiBulkDensity   , &    
    SoilResit4RootPentrate_vr       =>  plt_soilchem%SoilResit4RootPentrate_vr   , &    
    C4PhotosynDowreg_brch           =>  plt_photo%C4PhotosynDowreg_brch       &
  )  
  IF(SoiBulkDensity(L).GT.ZERO)THEN
    Root1stDepz2Surf=Root1stDepz_pft(N,NR,NZ)-CumSoilThickness(0)
  ELSE
    Root1stDepz2Surf=Root1stDepz_pft(N,NR,NZ)
  ENDIF

  IF(Root1stDepz2Surf.GT.CumSoilThickness(L-1).AND.ICHK1(N,NR).EQ.0)THEN
    Root1stXNumL_pvr(N,L,NZ)=Root1stXNumL_pvr(N,L,NZ)+RootAreaPopu
    IF(Root1stDepz2Surf.LE.CumSoilThickness(L).OR.L.EQ.MaxNumRootLays)THEN
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
!     SoilResit4RootPentrate_vr,SoilResit4PrimRootPentration=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
!     WFNR=water function for root extension
!     WFNRG=respiration function of root water potential
!
      SoilResit4PrimRootPentration=SoilResit4RootPentrate_vr(L)*Root1stRadius_pvr(N,L,NZ)/1.0E-03_r8
      WFNR=AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-PSIMin4OrganExtens-SoilResit4PrimRootPentration))
      IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
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
      IF(RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
        CNPG=AMIN1(RootNonstructElmConc_pvr(ielmn,N,L,NZ)/(RootNonstructElmConc_pvr(ielmn,N,L,NZ) &
          +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI),RootNonstructElmConc_pvr(ielmp,N,L,NZ) &
          /(RootNonstructElmConc_pvr(ielmp,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI))
      ELSE
        CNPG=1.0_r8
      ENDIF
!
!     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     Rmaint_CO2=root maintenance respiration
!     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT1N=primary root N mass
!     TFN6_vr=temperature function for root maintenance respiration
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     fRootGrowPSISense_vr=growth function of root water potential
!
      Rmaint_CO2=AZMAX1(RmSpecPlant*Root1stElm_raxs(ielmn,N,NR,NZ))*TFN6_vr(L)
      IF(is_root_shallow(iPlantRootProfile_pft(NZ)).OR.&
        iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
        Rmaint_CO2=Rmaint_CO2*fRootGrowPSISense_vr(N,L)
      ENDIF
!
!     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RNonstCO2_OUlm=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of primary root sink strength in axis
!     CPOOL=non-structural C mass
!     fTgrowRootP_vr=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!     fRootGrowPSISense_vr=growth function of root water potential
!
      RNonstCO2_OUlm=AZMAX1(VMXC*FRTN*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)) &
        *fTgrowRootP_vr(L,NZ)*CNPG*C4PhotosynDowreg_brch(MainBranchNum_pft(NZ),NZ)*fRootGrowPSISense_vr(N,L)
      IF(Root1stDepz2Surf.GE.CumSoilThickness(MaxNumRootLays))THEN
        RNonstCO2_OUlm=AMIN1(Rmaint_CO2,RNonstCO2_OUlm)
      ENDIF
!
!     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RAutoRootO2Limter_pvr=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RNonstCO2_Oltd=RNonstCO2_OUlm*RAutoRootO2Limter_pvr(N,L,NZ)
      RCO2XM=RNonstCO2_OUlm-Rmaint_CO2
      RCO2X=RNonstCO2_Oltd-Rmaint_CO2
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
!     RGrowCO2_OUlm,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))
      PPOOLB=AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ),PPOOLB*DMRTR/CPRTS(NZ))
      IF(RCO2YM.GT.0.0_r8)THEN
        RGrowCO2_OUlm=AMIN1(RCO2YM,FNP)
      ELSE
        RGrowCO2_OUlm=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0_r8)THEN
        RGrowCO2_Oltd=AMIN1(RCO2Y,FNP*RAutoRootO2Limter_pvr(N,L,NZ))
      ELSE
        RGrowCO2_Oltd=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     RootMycoNonst4Grow_Oulm,RootMycoNonst4Grow(ielmc)=total non-structural C used in growth and respn unltd,ltd by O2
!     RGrowCO2_OUlm,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,RootCYldO2ltd=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD1M,RootMycoNonst4Grow(ielmn),RootMycoNonst4Grow(ielmp)=nonstructural N,P unltd,ltd by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
      RootMycoNonst4Grow_Oulm=RGrowCO2_OUlm/DMRTD
      RootMycoNonst4Grow(ielmc)=RGrowCO2_Oltd/DMRTD
      GRTWGM=RootMycoNonst4Grow_Oulm*RootBiomGrowthYield(NZ)
      RootCYldO2ltd=RootMycoNonst4Grow(ielmc)*RootBiomGrowthYield(NZ)
      ZADD1M=AZMAX1(GRTWGM*CNRTW)
      RootMycoNonst4Grow(ielmn)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootCYldO2ltd*CNRTW))
      RootMycoNonst4Grow(ielmp)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootCYldO2ltd*CPRTW))
      CNRDM=AZMAX1(1.70_r8*ZADD1M)
      CNRDA=AZMAX1(1.70_r8*RootMycoNonst4Grow(ielmn))

      call PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     Rmaint_CO2=root maintenance respiration
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RootMycoNonst4Grow(ielmc)=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC1=fraction of primary root C to be remobilized
!     RootChemElmRemob(ielmc),RootChemElmRemob(ielmn),RootChemElmRemob(ielmp)=remobilization of C,N,P from senescing root
!     RootMycoNonst4Grow(ielmn),RootMycoNonst4Grow(ielmp)=nonstructural N,P ltd by O2 used in growth
!
      RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)&
          -AMIN1(Rmaint_CO2,RNonstCO2_Oltd)-CNRDA-SNCR

      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)&
          -RootMycoNonst4Grow(NE)+FSNC1*RootChemElmRemob(NE)
      ENDDO    
!
!     TOTAL PRIMARY ROOT RESPIRATION
!
!     RCO2T_Oulm,RCO2T_Oltd=total C respiration unltd,ltd by O2
!     Rmaint_CO2=root maintenance respiration
!     RNonstCO2_OUlm,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
!     RGrowCO2_OUlm,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A_pvr=total root respiration
!     RootRespPotent_pvr,RCO2N_pvr=RCO2A_pvr unltd by O2,nonstructural C
!
      RCO2T_Oulm=AMIN1(Rmaint_CO2,RNonstCO2_OUlm)+RGrowCO2_OUlm+SNCRM+CNRDM
      RCO2T_Oltd=AMIN1(Rmaint_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+SNCR+CNRDA
!
!     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
!     THROUGH WHICH PRIMARY ROOTS GROW
!
!     RTDP1=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     Root1stLen_rpvr=primary root length
!     SeedDepth_pft=seeding depth
!     FRCO2=fraction of primary root respiration attributed to layer
!     RCO2A_pvr=total root respiration
!     RootRespPotent_pvr,RCO2N_pvr=RCO2A_pvr unltd by O2,nonstructural C
!     RCO2T_Oulm,RCO2T_Oltd=total C respiration unltd,ltd by O2
!
      IF(Root1stDepz_pft(N,NR,NZ).GT.CumSoilThickness(NGTopRootLayer_pft(NZ)))THEN
        TFRCO2=0._r8
        D5100: DO LL=NGTopRootLayer_pft(NZ),NIXBotRootLayer_rpft(NR,NZ)
          IF(LL.LT.NIXBotRootLayer_rpft(NR,NZ))THEN
            FRCO2=AMIN1(1.0_r8,Root1stLen_rpvr(N,LL,NR,NZ)/(Root1stDepz_pft(N,NR,NZ)-SeedDepth_pft(NZ)))
          ELSE
            FRCO2=1.0_r8-TFRCO2
          ENDIF
          TFRCO2=TFRCO2+FRCO2
          RootRespPotent_pvr(N,LL,NZ)=RootRespPotent_pvr(N,LL,NZ)+RCO2T_Oulm*FRCO2
          RCO2N_pvr(N,LL,NZ)=RCO2N_pvr(N,LL,NZ)+RCO2T_Oltd*FRCO2
          RCO2A_pvr(N,LL,NZ)=RCO2A_pvr(N,LL,NZ)-RCO2T_Oltd*FRCO2
        ENDDO D5100
      ELSE
        RootRespPotent_pvr(N,L,NZ)=RootRespPotent_pvr(N,L,NZ)+RCO2T_Oulm
        RCO2N_pvr(N,L,NZ)=RCO2N_pvr(N,L,NZ)+RCO2T_Oltd
        RCO2A_pvr(N,L,NZ)=RCO2A_pvr(N,L,NZ)-RCO2T_Oltd
      ENDIF
!
!     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
!     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
!     HAVE DISAPPEARED
!
!     RootCYldO2ltd=primary root C growth ltd by O2
!     RootNetGrowthElms(ielmc),RootNetGrowthElms(ielmn),RootNetGrowthElms(ielmp)=net primary root C,N,P growth
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RootMycoNonst4Grow(ielmn),RootMycoNonst4Grow(ielmp)=nonstructural N,P ltd by O2 used in growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     Root2ndLen_pvr=secondary root length
!
      RootNetGrowthElms(ielmc)=RootCYldO2ltd-FSNC1*Root1stElm_raxs(ielmc,N,NR,NZ)
      RootNetGrowthElms(ielmn)=RootMycoNonst4Grow(ielmn)-FSNC1*Root1stElm_raxs(ielmn,N,NR,NZ)
      RootNetGrowthElms(ielmp)=RootMycoNonst4Grow(ielmp)-FSNC1*Root1stElm_raxs(ielmp,N,NR,NZ)
      IF(RootNetGrowthElms(ielmc).LT.0.0_r8)THEN
        LX=MAX(1,L-1)
        D5105: DO LL=L,LX,-1
          GRTWTM=RootNetGrowthElms(ielmc)
          DO NE=1,NumPlantChemElms
            IF(RootNetGrowthElms(NE).LT.0.0_r8)THEN
              IF(RootNetGrowthElms(NE).GT.-RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ))THEN
                if(NE==ielmc)Root2ndLen_pvr(N,LL,NR,NZ)=Root2ndLen_pvr(N,LL,NR,NZ)&
                  +RootNetGrowthElms(NE) &
                  *Root2ndLen_pvr(N,LL,NR,NZ)/RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)
                RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)&
                  +RootNetGrowthElms(NE)
                  RootNetGrowthElms(NE)=0._r8
              ELSE
                if(NE==ielmc)Root2ndLen_pvr(N,LL,NR,NZ)=0._r8
                RootNetGrowthElms(NE)=RootNetGrowthElms(NE)&
                  +RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)
                RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)=0._r8
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
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     WTRT2,WTRT2N,WTRT2P=mycorrhizal C,N,P mass
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     Root2ndLen_pvr=mycorrhizal length
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
!
          IF(GRTWTM.LT.0.0_r8)THEN
            IF(RootMyco2ndStrutElms_rpvr(ielmc,ipltroot,LL,NR,NZ).GT.ZEROP(NZ))THEN
              FSNCM=AMIN1(1.0_r8,ABS(GRTWTM)/RootMyco2ndStrutElms_rpvr(ielmc,ipltroot,LL,NR,NZ))
            ELSE
              FSNCM=1.0_r8
            ENDIF
            IF(RootMycoActiveBiomC_pvr(ipltroot,LL,NZ).GT.ZEROP(NZ))THEN
              FSNCP=AMIN1(1.0_r8,ABS(GRTWTM)/RootMycoActiveBiomC_pvr(ipltroot,LL,NZ))
            ELSE
              FSNCP=1.0_r8
            ENDIF

            
            D6450: DO M=1,jsken
              DO NE=1,NumPlantChemElms
                LitfalStrutElms_pvr(NE,M,k_woody_litr,LL,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,LL,NZ) &
                  +CFOPE(NE,icwood,M,NZ) &
                  *FSNCM*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ))*FWODRE(NE,k_woody_litr)

                LitfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ) &
                  +CFOPE(NE,iroot,M,NZ) &
                  *FSNCM*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ))*FWODRE(NE,k_fine_litr)

                LitfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ) &
                  +CFOPE(NE,inonstruct,M,NZ) &
                  *FSNCP*AZMAX1(RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ))
              ENDDO    
            ENDDO D6450
            DO NE=1,NumPlantChemElms  
              RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ)= &
                AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)

              RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ)=&
                AZMAX1(RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ))*(1.0_r8-FSNCP)

            ENDDO
            Root2ndLen_pvr(imycorrhz,LL,NR,NZ)=AZMAX1(Root2ndLen_pvr(imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)
          ENDIF
        ENDDO D5105
      ENDIF
!
      call PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,RootCYldO2ltd,RootNetGrowthElms,&
        TotRoot1stLen,Root1stC)
    ENDIF
!
!
    IF(L.EQ.NIXBotRootLayer_rpft(NR,NZ))THEN
      !bottom root layer
      call WithdrawPrimRoot(L,NR,NZ,N,Root1stDepz2Surf,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootAreaPopu)
    ENDIF

!
!     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
!
    IF(RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
        RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ) &
          +RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
      RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)=0._r8
    ENDIF
    IF(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
        RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ) &
          +RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
      RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)=0._r8
    ENDIF
!
!     TOTAL PRIMARY ROOT LENGTH AND MASS
!
!     TotRoot1stLen=total primary root length
!     Root1stC=total primary root C mass
!     Root1stLen_rpvr=primary root length in soil layer
!     WTRT1=primary root C mass in soil layer
!     NIXBotRootLayer_rpft=deepest root layer
!
    TotRoot1stLen=TotRoot1stLen+Root1stLen_rpvr(N,L,NR,NZ)
    Root1stC=Root1stC+RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
    NIXBotRootLayer_rpft(NR,NZ)=MIN(NIXBotRootLayer_rpft(NR,NZ),MaxNumRootLays)
    IF(L.EQ.NIXBotRootLayer_rpft(NR,NZ))NRX(N,NR)=1
  ENDIF
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
  real(r8) :: RootChemElmRemob(NumPlantChemElms)
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SNCR,SNCRM
! begin_execution
  associate(                          &
    Root1stElm_raxs                 =>  plt_biom%Root1stElm_raxs     , &
    RootNonstructElmConc_pvr        =>  plt_biom%RootNonstructElmConc_pvr     , &
    ZEROP                           =>  plt_biom%ZEROP      , &
    ZERO                            =>  plt_site%ZERO       , &
    icwood                          =>  pltpar%icwood       , &
    k_woody_litr                    => pltpar%k_woody_litr,&
    k_fine_litr                     => pltpar%k_fine_litr, &
    iroot                           =>  pltpar%iroot        , &
    FWODRE                          =>  plt_allom%FWODRE    , &
    CFOPE                           =>  plt_soilchem%CFOPE  , &
    RAutoRootO2Limter_pvr          =>  plt_rbgc%RAutoRootO2Limter_pvr       , &
    LitfalStrutElms_pvr               =>  plt_bgcr%LitfalStrutElms_pvr       , &
    iPlantCalendar_brch             =>  plt_pheno%iPlantCalendar_brch    , &
    MainBranchNum_pft               =>  plt_morph%MainBranchNum_pft         &
  )
!
!     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LitrFall
!
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 &
    .AND.RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmn,N,L,NZ)&
      ,RootNonstructElmConc_pvr(ielmn,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI) &
      ,safe_adb(RootNonstructElmConc_pvr(ielmp,N,L,NZ),RootNonstructElmConc_pvr(ielmp,N,L,NZ)&
        +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ),&
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ),&
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmp,N,L,NZ)/CPKI)))
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
!     RAutoRootO2Limter_pvr=constraint by O2 consumption on all root processes
!     RootChemElmRemob(ielmc),RootChemElmRemob(ielmn),RootChemElmRemob(ielmp)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC1=fraction of primary root C to be remobilized
!
  IF(-RCO2XM.GT.0.0_r8)THEN
    IF(-RCO2XM.LT.Root1stElm_raxs(ielmc,N,NR,NZ)*RCCC)THEN
      SNCRM=-RCO2XM
    ELSE
      SNCRM=AZMAX1(Root1stElm_raxs(ielmc,N,NR,NZ)*RCCC)
    ENDIF
  ELSE
    SNCRM=0._r8
  ENDIF
  IF(-RCO2X.GT.0.0_r8)THEN
    IF(-RCO2X.LT.Root1stElm_raxs(ielmc,N,NR,NZ)*RCCC)THEN
      SNCR=-RCO2X
    ELSE
      SNCR=AZMAX1(Root1stElm_raxs(ielmc,N,NR,NZ)*RCCC)*RAutoRootO2Limter_pvr(N,L,NZ)
    ENDIF
  ELSE
    SNCR=0._r8
  ENDIF
  IF(SNCR.GT.0.0_r8.AND.Root1stElm_raxs(ielmc,N,NR,NZ).GT.ZEROP(NZ))THEN
    RootChemElmRemob(ielmc)=RCCC*Root1stElm_raxs(ielmc,N,NR,NZ)
    RootChemElmRemob(ielmn)=Root1stElm_raxs(ielmn,N,NR,NZ)*(RCCN+(1.0_r8-RCCN) &
      *RootChemElmRemob(ielmc)/Root1stElm_raxs(ielmc,N,NR,NZ))
    RootChemElmRemob(ielmp)=Root1stElm_raxs(ielmp,N,NR,NZ)*(RCCP+(1.0_r8-RCCP) &
      *RootChemElmRemob(ielmc)/Root1stElm_raxs(ielmc,N,NR,NZ))
    IF(RootChemElmRemob(ielmc).GT.ZEROP(NZ))THEN
      FSNC1=AZMAX1(AMIN1(1.0_r8,SNCR/RootChemElmRemob(ielmc)))
    ELSE
      FSNC1=1.0_r8
    ENDIF
  ELSE
    RootChemElmRemob(1:NumPlantChemElms)=0._r8
    FSNC1=0._r8
  ENDIF
!
!     PRIMARY ROOT LitrFall CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RootChemElmRemob(ielmc),RootChemElmRemob(ielmn),RootChemElmRemob(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
  D6355: DO M=1,jsken
    DO NE=1,NumPlantChemElms    
      LitfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ) &
        +CFOPE(NE,icwood,M,NZ)*FSNC1*(Root1stElm_raxs(NE,N,NR,NZ)-RootChemElmRemob(NE))*FWODRE(NE,k_woody_litr)

      LitfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
        +CFOPE(NE,iroot,M,NZ)*FSNC1*(Root1stElm_raxs(NE,N,NR,NZ)-RootChemElmRemob(NE))*FWODRE(NE,k_fine_litr)
    ENDDO
  ENDDO D6355

  end associate
  end subroutine PrimRootRemobilization

!------------------------------------------------------------------------------------------
  subroutine PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,RootCYldO2ltd,RootNetGrowthElms,&
    TotRoot1stLen,Root1stC)
  implicit none
  integer, intent(in) :: L,L1
  integer, intent(in) :: N    !root/myco indicator
  integer, intent(in) :: NR   !root axis 
  integer, intent(in) :: NZ   !pft number
  real(r8), intent(in):: WFNR,FRTN
  real(r8), intent(in) :: RootCYldO2ltd !oxygen limited root C yield for growth
  real(r8), intent(in) :: RootNetGrowthElms(NumPlantChemElms)
  real(r8), intent(inout) :: TotRoot1stLen,Root1stC
  real(r8) :: Root1stExtension
  real(r8) :: FGROL,FGROZ
  integer :: NE
  REAL(R8) :: XFRE(NumPlantChemElms)
! begin_execution
  associate(                             &
    RootMyco1stStrutElms_rpvr        =>  plt_biom%RootMyco1stStrutElms_rpvr      , &
    Root1stElm_raxs              =>  plt_biom%Root1stElm_raxs      , &
    RootProteinC_pvr             =>  plt_biom%RootProteinC_pvr       , &
    PopuRootMycoC_pvr            =>  plt_biom% PopuRootMycoC_pvr      , &
    RootMycoNonstElms_rpvr         =>  plt_biom%RootMycoNonstElms_rpvr      , &
    ZEROP                        =>  plt_biom%ZEROP       , &
    rCNNonstructRemob_pft        =>  plt_allom%rCNNonstructRemob_pft      , &
    rCPNonstructRemob_pft        =>  plt_allom%rCPNonstructRemob_pft       , &
    FWODRE                       =>  plt_allom%FWODRE     , &
    CumSoilThickness             =>  plt_site%CumSoilThickness      , &
    DLYR3                        =>  plt_site%DLYR3       , &
    MaxNumRootLays               =>  plt_site%MaxNumRootLays          , &
    PlantPopulation_pft          =>  plt_site%PlantPopulation_pft          , &
    PSIRootTurg_vr               =>  plt_ew%PSIRootTurg_vr         , &
    PSIRootOSMO_vr               =>  plt_ew%PSIRootOSMO_vr         , &
    PSIRoot_pvr                  =>  plt_ew%PSIRoot_pvr        , &
    k_woody_litr                 =>  pltpar%k_woody_litr , &
    k_fine_litr                  =>  pltpar%k_fine_litr   , &
    Root1stSpecLen_pft           =>  plt_morph%Root1stSpecLen_pft     , &
    Root1stRadius_pvr            =>  plt_morph%Root1stRadius_pvr     , &
    Root1stLen_rpvr              =>  plt_morph%Root1stLen_rpvr      , &
    Root1stDepz_pft              =>  plt_morph%Root1stDepz_pft     , &
    NGTopRootLayer_pft           =>  plt_morph%NGTopRootLayer_pft        , &
    NIXBotRootLayer_rpft         =>  plt_morph%NIXBotRootLayer_rpft       , &
    SeedDepth_pft                =>  plt_morph%SeedDepth_pft        &
  )
!     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     Root1stExtension=primary root length extension
!     RootCYldO2ltd=primary root C growth ltd by O2
!     Root1stSpecLen_pft=specific primary root length from startq.f
!     PP=PFT population
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     RootNetGrowthElms(ielmc),RootNetGrowthElms(ielmn),RootNetGrowthElms(ielmp)=net primary root C,N,P growth
!     RTDP1=primary root depth from soil surface
!     SeedDepth_pft=seeding depth
!     FSNC1=fraction of primary root C to be remobilized
!     Root1stLen_rpvr=primary root length
!     RootNetGrowthElms(ielmc),RootNetGrowthElms(ielmn),RootNetGrowthElms(ielmp)=net root C,N,P growth
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     DLYR=soil layer thickness
!
  IF(RootNetGrowthElms(ielmc).LT.0.0_r8 .AND. Root1stElm_raxs(ielmc,N,NR,NZ).GT.ZEROP(NZ))THEN
    Root1stExtension=RootCYldO2ltd*Root1stSpecLen_pft(N,NZ)/PlantPopulation_pft(NZ)*WFNR*FWODRE(ielmc,k_fine_litr) &
      +RootNetGrowthElms(ielmc)*(Root1stDepz_pft(N,NR,NZ)-SeedDepth_pft(NZ))/Root1stElm_raxs(ielmc,N,NR,NZ)
  ELSE
    Root1stExtension=RootCYldO2ltd*Root1stSpecLen_pft(N,NZ)/PlantPopulation_pft(NZ)*WFNR*FWODRE(ielmc,k_fine_litr)
  ENDIF
  !the extension should not exceed soil layer thickness
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
  IF(Root1stExtension.GT.ZEROP(NZ) .AND. L.LT.MaxNumRootLays)THEN
    FGROL=AZMAX1(AMIN1(1.0_r8,(CumSoilThickness(L)-Root1stDepz_pft(N,NR,NZ))/Root1stExtension))
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
!     RootNetGrowthElms(ielmc),RootNetGrowthElms(ielmn),RootNetGrowthElms(ielmp)=net root C,N,P growth
!     Root1stExtension=primary root length extension
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     FGROL,FGROZ=fraction of Root1stExtension in current,next lower soil layer
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!     Root1stLen_rpvr=primary root length
!
  Root1stDepz_pft(N,NR,NZ)=Root1stDepz_pft(N,NR,NZ)+Root1stExtension

  DO NE=1,NumPlantChemElms
    Root1stElm_raxs(NE,N,NR,NZ)=Root1stElm_raxs(NE,N,NR,NZ)+RootNetGrowthElms(NE)
    RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
      +RootNetGrowthElms(NE)*FGROL
  ENDDO
  RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)+AMIN1(rCNNonstructRemob_pft(NZ)&
    *RootMyco1stStrutElms_rpvr(ielmn,N,L,NR,NZ) &
    ,rCPNonstructRemob_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmp,N,L,NR,NZ))
  Root1stLen_rpvr(N,L,NR,NZ)=Root1stLen_rpvr(N,L,NR,NZ)+Root1stExtension*FGROL
!
!     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
!     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
!     OF CURRENT SOIL LAYER
!
!     FGROZ=fraction of Root1stExtension in next lower soil layer
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     RootNetGrowthElms(ielmc),RootNetGrowthElms(ielmn),RootNetGrowthElms(ielmp)=net root C,N,P growth
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!     WTRTD=root C mass
!     Root1stLen_rpvr=primary root length
!     Root1stExtension=primary root length extension
!     FRTN=fraction of primary root sink strength in axis
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     PSIRoot_pvr,PSIRootTurg_vr,PSIRootOSMO_vr=root total,turgor,osmotic water potential
!     NIXBotRootLayer_rpft=deepest root layer
!
  IF(FGROZ.GT.0.0_r8)THEN
    DO NE=1,NumPlantChemElms
      RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ)&
        +RootNetGrowthElms(NE)*FGROZ
    ENDDO
    RootProteinC_pvr(N,L1,NZ)=RootProteinC_pvr(N,L1,NZ)+AMIN1(rCNNonstructRemob_pft(NZ)&
      *RootMyco1stStrutElms_rpvr(ielmn,N,L1,NR,NZ) &
      ,rCPNonstructRemob_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmp,N,L1,NR,NZ))
    PopuRootMycoC_pvr(N,L1,NZ)= PopuRootMycoC_pvr(N,L1,NZ)+RootMyco1stStrutElms_rpvr(ielmc,N,L1,NR,NZ)
    Root1stLen_rpvr(N,L1,NR,NZ)=Root1stLen_rpvr(N,L1,NR,NZ)+Root1stExtension*FGROZ
    Root1stRadius_pvr(N,L1,NZ)=Root1stRadius_pvr(N,L,NZ)
    TotRoot1stLen=TotRoot1stLen+Root1stLen_rpvr(N,L1,NR,NZ)
    Root1stC=Root1stC+RootMyco1stStrutElms_rpvr(ielmc,N,L1,NR,NZ)

    DO NE=1,NumPlantChemElms
      XFRE(NE)=FRTN* RootMycoNonstElms_rpvr(NE,N,L,NZ)
      RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)-XFRE(NE)
      RootMycoNonstElms_rpvr(NE,N,L1,NZ)=RootMycoNonstElms_rpvr(NE,N,L1,NZ)+XFRE(NE)
    ENDDO
    PSIRoot_pvr(N,L1,NZ)=PSIRoot_pvr(N,L,NZ)
    PSIRootOSMO_vr(N,L1,NZ)=PSIRootOSMO_vr(N,L,NZ)
    PSIRootTurg_vr(N,L1,NZ)=PSIRootTurg_vr(N,L,NZ)
    NIXBotRootLayer_rpft(NR,NZ)=MAX(NGTopRootLayer_pft(NZ),L+1)
  ENDIF
  end associate
  end subroutine PrimRootExtension
!------------------------------------------------------------------------------------------

  subroutine WithdrawPrimRoot(L,NR,NZ,N,Root1stDepz2Surf,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootAreaPopu)
  implicit none
  integer, intent(in) :: L,NR,NZ,N
  real(r8), intent(in):: Root1stDepz2Surf
  real(r8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: RootAreaPopu
  integer :: LL,NN,NE,NTG
  real(r8) :: XFRD,XFRW,FRTN
  real(r8) :: XFRE(NumPlantChemElms)
  logical :: RootDepzChk

! begin_execution
  associate(                             &
    RootMyco1stStrutElms_rpvr         =>  plt_biom%RootMyco1stStrutElms_rpvr   , &
    RootMyco2ndStrutElms_rpvr          =>  plt_biom%RootMyco2ndStrutElms_rpvr   , &
    RootMycoNonstElms_rpvr          =>  plt_biom%RootMycoNonstElms_rpvr   , &
    RootProteinC_pvr              =>  plt_biom%RootProteinC_pvr    , &
    PopuRootMycoC_pvr             =>  plt_biom% PopuRootMycoC_pvr   , &
    RootNodulStrutElms_pvr              =>  plt_biom%RootNodulStrutElms_pvr   , &
    ZEROP                         =>  plt_biom%ZEROP    , &
    RootNodulNonstElms_pvr         =>  plt_biom%RootNodulNonstElms_pvr  , &
    RootGasLossDisturb_pft        =>  plt_bgcr%RootGasLossDisturb_pft    , &
    CumSoilThickness              =>  plt_site%CumSoilThickness   , &
    ZEROS2                        =>  plt_site%ZEROS2   , &
    DLYR3                         =>  plt_site%DLYR3    , &
    trcg_rootml_pvr               =>  plt_rbgc%trcg_rootml_pvr , &
    trcs_rootml_pvr               =>  plt_rbgc%trcs_rootml_pvr, &
    VLSoilPoreMicP                =>  plt_soilchem%VLSoilPoreMicP , &
    Root1stXNumL_pvr              =>  plt_morph%Root1stXNumL_pvr    , &
    Root2ndXNum_pvr               =>  plt_morph%Root2ndXNum_pvr    , &
    Root1stLen_rpvr               =>  plt_morph%Root1stLen_rpvr   , &
    Root2ndXNum_rpvr              =>  plt_morph%Root2ndXNum_rpvr   , &
    Root1stDepz_pft               =>  plt_morph%Root1stDepz_pft  , &
    NGTopRootLayer_pft            =>  plt_morph%NGTopRootLayer_pft     , &
    iPlantNfixType                =>  plt_morph%iPlantNfixType  , &
    MY                            =>  plt_morph%MY      , &
    SeedDepth_pft                 =>  plt_morph%SeedDepth_pft   , &
    NIXBotRootLayer_rpft          =>  plt_morph%NIXBotRootLayer_rpft      &
  )
!     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
!     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
!     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
!     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
!
!     NIXBotRootLayer_rpft=deepest root layer
!     VLSoilPoreMicP=soil layer volume excluding macropore, rocks
!     Root1stDepz2Surf=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     SeedDepth_pft=seeding depth
!     FRTN=fraction of primary root sink strength in axis
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     Root1stLen_rpvr=primary root length
!     RootProteinC_pvr=root protein C mass
!     WTRTD=root C mass
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root

  D5115: DO LL=L,NGTopRootLayer_pft(NZ)+1,-1
    RootDepzChk=Root1stDepz2Surf.LT.CumSoilThickness(LL-1).OR.Root1stDepz2Surf.LT.SeedDepth_pft(NZ)
    IF(VLSoilPoreMicP(LL-1).GT.ZEROS2 .AND. RootDepzChk)THEN
      IF(RootSinkC_vr(N,LL).GT.ZEROP(NZ))THEN
        FRTN=(Root1stSink_pvr(N,LL,NR)+Root2ndSink_pvr(N,LL,NR))/RootSinkC_vr(N,LL)
      ELSE
        FRTN=1.0_r8
      ENDIF
      D5110: DO NN=1,MY(NZ)
        Root1stLen_rpvr(NN,LL-1,NR,NZ)=Root1stLen_rpvr(NN,LL-1,NR,NZ)+Root1stLen_rpvr(NN,LL,NR,NZ)
        Root1stLen_rpvr(NN,LL,NR,NZ)=0._r8
        DO NE=1,NumPlantChemElms
          RootMyco1stStrutElms_rpvr(NE,NN,LL-1,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,NN,LL-1,NR,NZ)&
            +RootMyco1stStrutElms_rpvr(NE,NN,LL,NR,NZ)
          RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ)&
            +RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ)
          RootMyco1stStrutElms_rpvr(NE,NN,LL,NR,NZ)=0._r8
          RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ)=0._r8
          XFRE(NE)=FRTN*RootMycoNonstElms_rpvr(NE,NN,LL,NZ)
          RootMycoNonstElms_rpvr(NE,NN,LL,NZ)=RootMycoNonstElms_rpvr(NE,NN,LL,NZ)-XFRE(NE)
          RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ)=RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ)+XFRE(NE)
        ENDDO
        XFRW=FRTN*RootProteinC_pvr(NN,L,NZ)
        XFRD=FRTN* PopuRootMycoC_pvr(NN,LL,NZ)

        RootProteinC_pvr(NN,LL,NZ)=RootProteinC_pvr(NN,LL,NZ)-XFRW
        RootProteinC_pvr(NN,LL-1,NZ)=RootProteinC_pvr(NN,LL-1,NZ)+XFRW
        PopuRootMycoC_pvr(NN,LL,NZ)= PopuRootMycoC_pvr(NN,LL,NZ)-XFRD
        PopuRootMycoC_pvr(NN,LL-1,NZ)= PopuRootMycoC_pvr(NN,LL-1,NZ)+XFRD
!
!     WITHDRAW GASES IN PRIMARY ROOTS
!
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
!     FRTN=fraction of primary root sink strength in axis
!
        DO NTG=idg_beg,idg_end-1
          RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-FRTN &
            *(trcg_rootml_pvr(idg_CO2,NN,LL,NZ)+trcs_rootml_pvr(idg_CO2,NN,LL,NZ))
          trcg_rootml_pvr(NTG,NN,LL,NZ)=(1.0_r8-FRTN)*trcg_rootml_pvr(NTG,NN,LL,NZ)
          trcs_rootml_pvr(NTG,NN,LL,NZ)=(1.0_r8-FRTN)*trcs_rootml_pvr(NTG,NN,LL,NZ)
        ENDDO
      ENDDO D5110
!
!     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
!
!     Root2ndXNum_rpvr,Root2ndXNum_pvr=number of secondary root axes
!     Root1stXNumL_pvr=number of primary root axes
!     Root1stLen_rpvr=primary root length
!     CumSoilThickness=depth from soil surface to layer bottom
!     SeedDepth_pft=seeding depth
!
      Root2ndXNum_pvr(N,LL,NZ)=Root2ndXNum_pvr(N,LL,NZ)-Root2ndXNum_rpvr(N,LL,NR,NZ)
      Root2ndXNum_pvr(N,LL-1,NZ)=Root2ndXNum_pvr(N,LL-1,NZ)+Root2ndXNum_rpvr(N,LL,NR,NZ)
      Root2ndXNum_rpvr(N,LL,NR,NZ)=0._r8
      Root1stXNumL_pvr(N,LL,NZ)=Root1stXNumL_pvr(N,LL,NZ)-RootAreaPopu
      IF(LL-1.GT.NGTopRootLayer_pft(NZ))THEN
        Root1stLen_rpvr(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness(LL-1)-Root1stDepz_pft(N,NR,NZ))
      ELSE
        Root1stLen_rpvr(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness(LL-1)-Root1stDepz_pft(N,NR,NZ)) &
          -(SeedDepth_pft(NZ)-CumSoilThickness(LL-2))
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
        DO NE=1,NumPlantChemElms
          XFRE(NE)=FRTN*RootNodulStrutElms_pvr(NE,LL,NZ)
          RootNodulStrutElms_pvr(NE,LL,NZ)=RootNodulStrutElms_pvr(NE,LL,NZ)-XFRE(NE)
          RootNodulStrutElms_pvr(NE,LL-1,NZ)=RootNodulStrutElms_pvr(NE,LL-1,NZ)+XFRE(NE)
          XFRE(NE)=FRTN*RootNodulNonstElms_pvr(NE,LL,NZ)
          RootNodulNonstElms_pvr(NE,LL,NZ)=RootNodulNonstElms_pvr(NE,LL,NZ)-XFRE(NE)
          RootNodulNonstElms_pvr(NE,LL-1,NZ)=RootNodulNonstElms_pvr(NE,LL-1,NZ)+XFRE(NE)
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

  subroutine NonstructlBiomTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,&
    Root2ndSink_pvr,RootSinkC,BegRemoblize)
  !
  !DESCRIPTION
  !transfer of nonstructural C/N/P 
  implicit none
  integer,  intent(in) :: I,J,NZ,BegRemoblize
  real(r8), intent(in):: PTRT
  real(r8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: RootSinkC(jroots)
  integer :: L,NB,N,NR,NE
  real(r8) :: ZPOOLS,ZPOOLT
  real(r8) :: ZPOOLB
  real(r8) :: ZPOOLD,EPOOLD
  real(r8) :: XFREX(NumPlantChemElms)
  real(r8) :: TotNonstElm_loc(NumPlantChemElms)
  real(r8) :: WTLSBZ(NumOfCanopyLayers1)
  real(r8) :: NonstElm_loc(NumPlantChemElms,NumOfCanopyLayers1)
  REAL(R8) :: FWTR(JZ1),FWTB(JP1)
  real(r8) :: CPOOLT
  real(r8) :: NonstElmRootE,NonstElmBrchE
  real(r8) :: CNL,CPL,NonstElmGradt
  real(r8) :: CPOOLB,CPOOLS
  real(r8) :: FWTC
  real(r8) :: FWTS
  real(r8) :: PPOOLB
  real(r8) :: PPOOLT
  real(r8) :: PTSHTR
  real(r8) :: PPOOLS
  real(r8) :: WTPLTT
  real(r8) :: WTRTLX
  real(r8) :: WTRTD1
  real(r8) :: WTSTKT
  real(r8) :: TotStalkRsrv_loc(NumPlantChemElms)
  real(r8) :: StalkRsrvGradt
  real(r8) :: WTRTD2,WTLSBX,WTLSBB
  real(r8) :: WTRTLR
  real(r8) :: XFRE(NumPlantChemElms)
  real(r8) :: sumchk1,sumchk2
!     begin_execution
  associate(                                &
    FWODBE                             =>   plt_allom%FWODBE  , &
    FWODRE                             =>   plt_allom%FWODRE  , &
    RootNonstructElmConc_pvr           =>   plt_biom%RootNonstructElmConc_pvr   , &
    Root1stElm_raxs                    =>   plt_biom%Root1stElm_raxs   , &
    RootMyco1stStrutElms_rpvr              =>   plt_biom%RootMyco1stStrutElms_rpvr   , &
    RootMyco2ndStrutElms_rpvr               =>   plt_biom%RootMyco2ndStrutElms_rpvr   , &
    RootMycoActiveBiomC_pvr            =>   plt_biom%RootMycoActiveBiomC_pvr   , &
    RootMycoNonstElms_rpvr               =>   plt_biom%RootMycoNonstElms_rpvr   , &
    LeafPetolBiomassC_brch             =>   plt_biom%LeafPetolBiomassC_brch   , &
    CanopyNonstElms_brch                  =>   plt_biom%CanopyNonstElms_brch   , &
    PopuRootMycoC_pvr                  =>   plt_biom% PopuRootMycoC_pvr   , &
    StalkBiomassC_brch                 =>   plt_biom%StalkBiomassC_brch   , &
    StalkRsrvElms_brch                 =>   plt_biom%StalkRsrvElms_brch  , &
    SeasonalNonstElms_pft                  =>   plt_biom%SeasonalNonstElms_pft    , &
    CanopyLeafShethC_pft               =>   plt_biom%CanopyLeafShethC_pft     , &
    RootElms_pft                       =>   plt_biom%RootElms_pft    , &
    ZEROL                              =>   plt_biom%ZEROL    , &
    ZEROP                              =>   plt_biom%ZEROP    , &
    iPlantTurnoverPattern_pft          =>   plt_pheno%iPlantTurnoverPattern_pft  , &
    iPlantRootProfile_pft              =>   plt_pheno%iPlantRootProfile_pft  , &
    iPlantBranchState_brch             =>   plt_pheno%iPlantBranchState_brch   , &
    iPlantPhenolPattern_pft            =>   plt_pheno%iPlantPhenolPattern_pft  , &
    ShutRutNonstructElmntConducts_pft  =>   plt_pheno%ShutRutNonstructElmntConducts_pft  , &
    iPlantCalendar_brch                =>   plt_pheno%iPlantCalendar_brch  , &
    Hours2LeafOut_brch                 =>   plt_pheno%Hours2LeafOut_brch    , &
    RCO2A_pvr                          =>   plt_rbgc%RCO2A_pvr    , &
    GrossResp_pft                      =>   plt_bgcr%GrossResp_pft    , &
    ECO_ER_col                         =>   plt_bgcr%ECO_ER_col     , &
    Eco_AutoR_col                      =>   plt_bgcr%Eco_AutoR_col     , &
    NU                                 =>   plt_site%NU       , &
    ZERO                               =>   plt_site%ZERO     , &
    k_woody_litr                       =>   pltpar%k_woody_litr,&
    k_fine_litr                        =>   pltpar%k_fine_litr, &
    NIXBotRootLayer_pft                =>   plt_morph%NIXBotRootLayer_pft    , &
    NIXBotRootLayer_rpft               =>   plt_morph%NIXBotRootLayer_rpft    , &
    Radius2ndRoot_pvr                  =>   plt_morph%Radius2ndRoot_pvr   , &
    MaxSoiL4Root                       =>   plt_morph%MaxSoiL4Root      , &
    MY                                 =>   plt_morph%MY      , &
    NumRootAxes_pft                    =>   plt_morph%NumRootAxes_pft    , &
    NumOfBranches_pft                  =>   plt_morph%NumOfBranches_pft       &
  )
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH LEAVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!     WHEN SEASONAL STORAGE C IS NOT BEING MOBILIZED
!
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     Hours2LeafOut_brch=hourly leafout counter
!     HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
!     LeafPetolBiomassC_brch=leaf+petiole mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!
  IF(NumOfBranches_pft(NZ).GT.1)THEN
    WTPLTT=0._r8
    TotNonstElm_loc(1:NumPlantChemElms)=0._r8
    D300: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
        IF(Hours2LeafOut_brch(NB,NZ).GT.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)))THEN
          WTLSBZ(NB)=AZMAX1(LeafPetolBiomassC_brch(NB,NZ))
          DO NE=1,NumPlantChemElms
            NonstElm_loc(NE,NB)=AZMAX1(CanopyNonstElms_brch(NE,NB,NZ))
            TotNonstElm_loc(ielmc)=TotNonstElm_loc(ielmc)+NonstElm_loc(ielmc,NB)
          ENDDO
          WTPLTT=WTPLTT+WTLSBZ(NB)
        ENDIF
      ENDIF
    ENDDO D300

    !!Nonst check

    !
    D305: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
        !leaf out criterion met
        IF(Hours2LeafOut_brch(NB,NZ).GT.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)))THEN
          IF(WTPLTT.GT.ZEROP(NZ) .AND. TotNonstElm_loc(ielmc).GT.ZEROP(NZ))THEN
            NonstElmGradt=TotNonstElm_loc(ielmc)*WTLSBZ(NB)-NonstElm_loc(ielmc,NB)*WTPLTT
            XFRE(ielmc)=0.01_r8*NonstElmGradt/WTPLTT
            DO NE=2,NumPlantChemElms
              NonstElmGradt=TotNonstElm_loc(NE)*NonstElm_loc(ielmc,NB)-NonstElm_loc(NE,NB)*TotNonstElm_loc(ielmc)
              XFRE(NE)=0.01_r8*NonstElmGradt/TotNonstElm_loc(ielmc)
            ENDDO
            DO NE=1,NumPlantChemElms
              CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
              if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8)then
                write(*,*)'1750CanopyNonstElms_brch(NE,NB,NZ)',NE,NB,CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE),XFRE(NE)
                stop
              endif
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
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     StalkBiomassC_brch=stalk sapwood mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
! the algorithm below is different from that for within canopy transfer between different branches
! why?
  IF(NumOfBranches_pft(NZ).GT.1)THEN
    WTSTKT=0._r8
    TotStalkRsrv_loc(1:NumPlantChemElms)=0._r8
    D330: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
        IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
          WTSTKT=WTSTKT+StalkBiomassC_brch(NB,NZ)
          DO NE=1,NumPlantChemElms
            TotStalkRsrv_loc(NE)=TotStalkRsrv_loc(NE)+StalkRsrvElms_brch(NE,NB,NZ)
          ENDDO
        ENDIF
      ENDIF
    ENDDO D330
    sumchk1=TotStalkRsrv_loc(ielmc)
    sumchk2=0._r8
    IF(WTSTKT.GT.ZEROP(NZ).AND.TotStalkRsrv_loc(ielmc).GT.ZEROP(NZ))THEN
      D335: DO NB=1,NumOfBranches_pft(NZ)
        IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
          IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
            StalkRsrvGradt=TotStalkRsrv_loc(ielmc)*StalkBiomassC_brch(NB,NZ)-StalkRsrvElms_brch(ielmc,NB,NZ)*WTSTKT
            XFRE(ielmc)=0.1_r8*StalkRsrvGradt/WTSTKT            
            StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
            sumchk2=sumchk2+StalkRsrvElms_brch(ielmc,NB,NZ)
            !based on stoichiometry gradient
            DO NE=2,NumPlantChemElms
              StalkRsrvGradt=TotStalkRsrv_loc(NE)*StalkRsrvElms_brch(ielmc,NB,NZ) &
                -StalkRsrvElms_brch(NE,NB,NZ)*TotStalkRsrv_loc(ielmc)
              XFRE(NE)=0.1_r8*StalkRsrvGradt/TotStalkRsrv_loc(ielmc)
              StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
            ENDDO
          ENDIF
        ENDIF
      ENDDO D335
    ENDIF
    print*,'1818sumchk stalkrsrv',sumchk1,sumchk2
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
  !this enables extension to other mycorrhizae
  IF(MY(NZ).GT.ipltroot)THEN
    DO N=2,MY(NZ)
      D425: DO L=NU,NIXBotRootLayer_pft(NZ)
        IF(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ).GT.ZEROP(NZ).AND. PopuRootMycoC_pvr(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
          !root
          WTRTD1=PopuRootMycoC_pvr(ipltroot,L,NZ)
          WTRTD2=AMIN1(PopuRootMycoC_pvr(ipltroot,L,NZ),AMAX1(FSNK &
            *PopuRootMycoC_pvr(ipltroot,L,NZ), PopuRootMycoC_pvr(N,L,NZ)))
          WTPLTT=WTRTD1+WTRTD2
          IF(WTPLTT.GT.ZEROP(NZ))THEN
            NonstElmGradt=(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*WTRTD2 &
              - RootMycoNonstElms_rpvr(ielmc,N,L,NZ)*WTRTD1)/WTPLTT
            XFRE(ielmc)=FMYC*NonstElmGradt
            RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
            RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+XFRE(ielmc)

            CPOOLT=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
            IF(CPOOLT.GT.ZEROP(NZ))THEN
              !exchange based on stoichiometry gradient
              DO NE=2,NumPlantChemElms
                NonstElmGradt=(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ) &
                  - RootMycoNonstElms_rpvr(NE,N,L,NZ)*RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/CPOOLT
                XFRE(NE)=FMYC*NonstElmGradt
                RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE(NE)
                RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+XFRE(NE)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO D425
    ENDDO
  ENDIF
!
!     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
!     IN PERENNIALS
!
  IF(BegRemoblize.EQ.itrue .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
    D5545: DO N=1,MY(NZ)
      D5550: DO L=NU,MaxSoiL4Root(NZ)
        IF(RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
          CNL=RootNonstructElmConc_pvr(ielmc,N,L,NZ)/(RootNonstructElmConc_pvr(ielmc,N,L,NZ) &
            +RootNonstructElmConc_pvr(ielmn,N,L,NZ)/CNKI)
          CPL=RootNonstructElmConc_pvr(ielmc,N,L,NZ)/(RootNonstructElmConc_pvr(ielmc,N,L,NZ) &
            +RootNonstructElmConc_pvr(ielmp,N,L,NZ)/CPKI)
        ELSE
          CNL=0._r8
          CPL=0._r8
        ENDIF
        XFREX(ielmc)=FXFR(iPlantTurnoverPattern_pft(NZ))*AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
        XFREX(ielmn)=FXFR(iPlantTurnoverPattern_pft(NZ))*AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*(1.0_r8+CNL)
        XFREX(ielmp)=FXFR(iPlantTurnoverPattern_pft(NZ))*AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*(1.0_r8+CPL)
        XFRE(ielmc)=AMIN1(XFREX(ielmc),XFREX(ielmn)/CNMN,XFREX(ielmp)/CPMN)
        XFRE(ielmn)=AMIN1(XFREX(ielmn),XFRE(ielmc)*CNMX,XFREX(ielmp)*CNMX/CPMN*0.5_r8)
        XFRE(ielmp)=AMIN1(XFREX(ielmp),XFRE(ielmc)*CPMX,XFREX(ielmn)*CPMX/CNMN*0.5_r8)
        DO NE=1,NumPlantChemElms
          RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)-XFRE(NE)
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+XFRE(NE)
          if(RootMycoNonstElms_rpvr(NE,N,L,NZ)<0._r8 .or. SeasonalNonstElms_pft(NE,NZ)<0._r8)then
            write(*,*)'1871RootMycoNonstElms_rpvr(NE,N,L,NZ)',NE,N,RootMycoNonstElms_rpvr(NE,N,L,NZ)+XFRE(NE),&
              SeasonalNonstElms_pft(NE,NZ)-XFRE(NE),XFRE(NE)
            stop
          endif
        ENDDO

      ENDDO D5550
    ENDDO D5545
  ENDIF
!
!     ROOT AND NODULE TOTALS
!
!     RootMycoActiveBiomC_pvr,WTRTD=active,actual root C mass
!     WTRT1,WTRT2=primary,secondary root C mass in soil layer
!     GrossResp_pft=total PFT respiration
!     RCO2A_pvr=total root respiration
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!
  D5445: DO N=1,MY(NZ)
    D5450: DO L=NU,MaxSoiL4Root(NZ)
      RootMycoActiveBiomC_pvr(N,L,NZ)=0._r8
      PopuRootMycoC_pvr(N,L,NZ)=0._r8
      D5460: DO NR=1,NumRootAxes_pft(NZ)
        RootMycoActiveBiomC_pvr(N,L,NZ)=RootMycoActiveBiomC_pvr(N,L,NZ)+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
        PopuRootMycoC_pvr(N,L,NZ)=PopuRootMycoC_pvr(N,L,NZ)+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ) &
          +RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
      ENDDO D5460
      GrossResp_pft(NZ)=GrossResp_pft(NZ)+RCO2A_pvr(N,L,NZ)
      ECO_ER_col=ECO_ER_col+RCO2A_pvr(N,L,NZ)
      Eco_AutoR_col=Eco_AutoR_col+RCO2A_pvr(N,L,NZ)
    ENDDO D5450

    DO  NR=1,NumRootAxes_pft(NZ)
      RootMycoActiveBiomC_pvr(N,NIXBotRootLayer_rpft(NR,NZ),NZ)=RootMycoActiveBiomC_pvr(N,NIXBotRootLayer_rpft(NR,NZ),NZ)&
        +Root1stElm_raxs(ielmc,N,NR,NZ)
    ENDDO
  ENDDO D5445
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
!
!     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
!     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     WTLS,WTRT=total PFT leaf+petiole,root C mass
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     RootSinkC_vr,RootSinkC=root layer,root system sink strength
!
!     IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
!  why 2/3 here?
  IF(CanopyLeafShethC_pft(NZ).GT.ZEROP(NZ))THEN
    FWTC=AMIN1(1.0_r8,0.667_r8*RootElms_pft(ielmc,NZ)/CanopyLeafShethC_pft(NZ))
  ELSE
    FWTC=1.0_r8
  ENDIF
  IF(RootElms_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
    FWTS=AMIN1(1.0_r8,CanopyLeafShethC_pft(NZ)/(0.667_r8*RootElms_pft(ielmc,NZ)))
  ELSE
    FWTS=1.0_r8
  ENDIF
!     ELSE
!     FWTC=1.0_r8
!     FWTS=1.0_r8
!     ENDIF
  D290: DO L=NU,MaxSoiL4Root(NZ)
    IF(RootSinkC(ipltroot).GT.ZEROP(NZ))THEN
      FWTR(L)=AZMAX1(RootSinkC_vr(ipltroot,L)/RootSinkC(ipltroot))
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
    CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)
  ENDDO D309
  
!  IF(NZ==1)THEN
!  write(153,*)I+J/24.,CanopyLeafShethC_pft(NZ)
!  ELSE
!  write(154,*)I+J/24.,CanopyLeafShethC_pft(NZ)
!  ENDIF
  
!
!     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
!     OF TOTAL SINK STRENGTH OF THE CANOPY
!
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date for setting final seed number
!     FWTB=branch sink weighting factor
!     ShutRutNonstructElmntConducts_pft=rate constant for equilibrating shoot-root nonstructural C concn from PFT file
!     PTRT=allocation to leaf+petiole used to modify ShutRutNonstructElmntConducts_pftin annuals
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     FWODB=C woody fraction in branch:0=woody,1=non-woody
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  D310: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      IF(CanopyLeafShethC_pft(NZ).GT.ZEROP(NZ))THEN
        FWTB(NB)=AZMAX1(LeafPetolBiomassC_brch(NB,NZ)/CanopyLeafShethC_pft(NZ))
      ELSE
        FWTB(NB)=1.0_r8
      ENDIF
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        PTSHTR=ShutRutNonstructElmntConducts_pft(NZ)*PTRT**0.167_r8
      ELSE
        PTSHTR=ShutRutNonstructElmntConducts_pft(NZ)
      ENDIF

      D415: DO L=NU,MaxSoiL4Root(NZ)
        WTLSBX=LeafPetolBiomassC_brch(NB,NZ)*FWODBE(ielmc,k_fine_litr)*FWTR(L)*FWTC
        WTRTLX=RootMycoActiveBiomC_pvr(ipltroot,L,NZ)*FWODRE(ielmc,k_fine_litr)*FWTB(NB)*FWTS
        WTLSBB=AZMAX1(WTLSBX,FSNK*WTRTLX)
        WTRTLR=AZMAX1(WTRTLX,FSNK*WTLSBX)
        WTPLTT=WTLSBB+WTRTLR
        IF(WTPLTT.GT.ZEROP(NZ))THEN
          CPOOLB=AZMAX1(CanopyNonstElms_brch(ielmc,NB,NZ)*FWTR(L))
          CPOOLS=AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*FWTB(NB))
          NonstElmGradt=(CPOOLB*WTRTLR-CPOOLS*WTLSBB)/WTPLTT
          XFRE(ielmc)=PTSHTR*NonstElmGradt
          CPOOLT=CPOOLS+CPOOLB
          CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)-XFRE(ielmc)
          RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+XFRE(ielmc)

          !N & P tranfer based on stoichiometry ratio
          IF(CPOOLT.GT.ZEROP(NZ))THEN
            DO NE=2,NumPlantChemElms
              NonstElmBrchE=AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ)*FWTR(L))
              NonstElmRootE=AZMAX1(RootMycoNonstElms_rpvr(ielmn,ipltroot,L,NZ)*FWTB(NB))
              NonstElmGradt=(NonstElmBrchE*CPOOLS-NonstElmRootE*CPOOLB)/CPOOLT
              XFRE(NE)=PTSHTR*NonstElmGradt
              IF(XFRE(NE)>0._r8)then
                XFRE(NE)=AMIN1(CanopyNonstElms_brch(NE,NB,NZ),XFRE(NE))
              ELSE  
                XFRE(NE)=AMAX1(-RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ),XFRE(NE))
              ENDIF
              CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
              RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)+XFRE(NE)              
              if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8 .or. RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)<0._r8)then
                write(*,*)'2012CanopyNonstElms_brch(NE,NB,NZ)',NE,NB,NZ,CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)&
                  ,RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE(NE) &
                  ,XFRE(NE),NonstElmBrchE,NonstElmRootE,NonstElmGradt,FWTR(L),FWTB(NB),PTSHTR
                stop
              endif
            ENDDO
          ENDIF  
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
  REAL(R8) :: Root1stLocDepz_vr(NumOfCanopyLayers1,JZ1)
  real(r8) :: RecoRootMycoC4Nup,CUPRO,CUPRC
  real(r8) :: RTDPP,RTDPS,RTSKP
  real(r8) :: RTSKS

  associate(                              &
    RootMycoNonstElms_rpvr      =>   plt_biom%RootMycoNonstElms_rpvr    , &
    ZEROP                     =>   plt_biom%ZEROP     , &
    iPlantRootProfile_pft     =>   plt_pheno%iPlantRootProfile_pft   , &
    VLSoilPoreMicP            =>   plt_soilchem%VLSoilPoreMicP  , &
    ZEROS2                    =>   plt_site%ZEROS2    , &
    NU                        =>   plt_site%NU        , &
    CumSoilThickness          =>   plt_site%CumSoilThickness    , &
    ZERO                      =>   plt_site%ZERO      , &
    DLYR3                     =>   plt_site%DLYR3     , &
    RootNutUptake_pvr         =>   plt_rbgc%RootNutUptake_pvr    , &
    RootOUlmNutUptake_pvr     =>   plt_rbgc%RootOUlmNutUptake_pvr    , &
    RootCUlmNutUptake_pvr     =>   plt_rbgc%RootCUlmNutUptake_pvr    , &
    RootMycoExudElm_pvr       =>   plt_rbgc%RootMycoExudElm_pvr    , &
    RCO2N_pvr                 =>   plt_rbgc%RCO2N_pvr     , &
    RootRespPotent_pvr        =>   plt_rbgc%RootRespPotent_pvr     , &
    RCO2A_pvr                 =>   plt_rbgc%RCO2A_pvr     , &
    CanPHeight4WatUptake      =>   plt_morph%CanPHeight4WatUptake    , &
    MY                        =>   plt_morph%MY       , &
    Root1stRadius_pvr         =>   plt_morph%Root1stRadius_pvr   , &
    Root1stDepz_pft           =>   plt_morph%Root1stDepz_pft   , &
    HypoctoHeight_pft         =>   plt_morph%HypoctoHeight_pft    , &
    Radius2ndRoot_pvr         =>   plt_morph%Radius2ndRoot_pvr    , &
    Root2ndXNum_rpvr          =>   plt_morph%Root2ndXNum_rpvr    , &
    AveLen2ndRoot_pvr         =>   plt_morph%AveLen2ndRoot_pvr    , &
    SeedDepth_pft             =>   plt_morph%SeedDepth_pft    , &
    MaxSoiL4Root              =>   plt_morph%MaxSoiL4Root       , &
    NumRootAxes_pft           =>   plt_morph%NumRootAxes_pft       &
  )

!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER

  RootSinkC_vr=0._R8
  Root1stSink_pvr=0._r8
  Root2ndSink_pvr=0._r8
  RootSinkC=0._r8
  D4995: DO N=1,MY(NZ)
    D4990: DO L=NU,MaxSoiL4Root(NZ)
!
!     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
!     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
!
!     VLSoilPoreMicP=soil layer volume excluding macropore, rocks
!     RecoRootMycoC4Nup=C respiration for nutrient uptake
!     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
!     RootNH4Uptake_pvr,RootNH4BUptake_pvr,RUPN03,RootNO3BUptake_pvr=uptake from non-band,band of NH4,NO3
!     RootH2PO4Uptake_pvr,RootH2PO4BUptake_pvr,RootHPO4Uptake_pvr,RootNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4
!     RootOUlmNutUptake_pvr,RootOUlmNutUptake_pvr,RUON03,RootOUlmNutUptake_pvr=uptake from non-band,band of NH4,NO3 unlimited by O2
!     RootOUlmNutUptake_pvr,RootOUlmNutUptake_pvr,RootOUlmNutUptake_pvr(ids_H1PO4,,RootOUlmNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4 unlimited by O2
!     RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr,RUCN03,RootCUlmNutUptake_pvr=uptake from non-band,band of NH4,NO3 unlimited by nonstructural C
!     RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4 unlimited by nonstructural C
!     why is 0.86? it refers to C cost for N asimilation
      IF(VLSoilPoreMicP(L).GT.ZEROS2)THEN
        RecoRootMycoC4Nup=0.86_r8*(RootNutUptake_pvr(ids_NH4,N,L,NZ)+RootNutUptake_pvr(ids_NH4B,N,L,NZ) &
          +RootNutUptake_pvr(ids_NO3,N,L,NZ)+RootNutUptake_pvr(ids_NO3B,N,L,NZ) &
          +RootNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootNutUptake_pvr(ids_H2PO4B,N,L,NZ) &
          +RootNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ))
        CUPRO=0.86_r8*(RootOUlmNutUptake_pvr(ids_NH4,N,L,NZ)+RootOUlmNutUptake_pvr(ids_NH4B,N,L,NZ) &
          +RootOUlmNutUptake_pvr(ids_NO3,N,L,NZ)+RootOUlmNutUptake_pvr(ids_NO3B,N,L,NZ)&
          +RootOUlmNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootOUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ)&
          +RootOUlmNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootOUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ))
        CUPRC=0.86_r8*(RootCUlmNutUptake_pvr(ids_NH4,N,L,NZ)+RootCUlmNutUptake_pvr(ids_NH4B,N,L,NZ) &
          +RootCUlmNutUptake_pvr(ids_NO3,N,L,NZ)+RootCUlmNutUptake_pvr(ids_NO3B,N,L,NZ) &
          +RootCUlmNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootCUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ)&
          +RootCUlmNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootCUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ))
!
!     ACCUMULATE RESPIRATION IN FLUX ARRAYS
!
!     RCO2A_pvr=total root respiration
!     RootRespPotent_pvr,RCO2N_pvr=RCO2A_pvr unltd by O2,nonstructural C
!     RecoRootMycoC4Nup=C respiration for nutrient uptake
!     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
!     CPOOLR=non-structural C mass in root
!
        RootRespPotent_pvr(N,L,NZ)=RootRespPotent_pvr(N,L,NZ)+CUPRO
        RCO2N_pvr(N,L,NZ)=RCO2N_pvr(N,L,NZ)+CUPRC
        RCO2A_pvr(N,L,NZ)=RCO2A_pvr(N,L,NZ)-RecoRootMycoC4Nup
        RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RecoRootMycoC4Nup
        print*,'2136RecoRootMycoC4Nup',RecoRootMycoC4Nup
!
!     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
!     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
!     RootNH4Uptake_pvr,RootNH4BUptake_pvr,RUPN03,RootNO3BUptake_pvr=uptake from non-band,band of NH4,NO3
!     RootH2PO4Uptake_pvr,RootH2PO4BUptake_pvr,RootHPO4Uptake_pvr,RootNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4
!
        D195: DO K=1,jcplx
          DO NE=1,NumPlantChemElms
             RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMycoExudElm_pvr(NE,N,K,L,NZ)
          ENDDO
        ENDDO D195
 
        if(RootMycoNonstElms_rpvr(ielmn,N,L,NZ)<0._r8 .or. RootMycoNonstElms_rpvr(ielmp,N,L,NZ)<0._r8)stop 
        RootMycoNonstElms_rpvr(ielmn,N,L,NZ)=RootMycoNonstElms_rpvr(ielmn,N,L,NZ)+&
          (RootNutUptake_pvr(ids_NH4,N,L,NZ) &
           +RootNutUptake_pvr(ids_NH4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_NO3,N,L,NZ)+RootNutUptake_pvr(ids_NO3B,N,L,NZ))
         RootMycoNonstElms_rpvr(ielmp,N,L,NZ)=RootMycoNonstElms_rpvr(ielmp,N,L,NZ) &
           +(RootNutUptake_pvr(ids_H2PO4,N,L,NZ)&
           +RootNutUptake_pvr(ids_H2PO4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ))
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
!     RRAD1,Radius2ndRoot_pvr=primary, secondary root radius
!     RootSinkC,RootSinkC_vr=total root sink strength
!
          IF(N.EQ.ipltroot)THEN
            IF(Root1stDepz_pft(N,NR,NZ).GT.CumSoilThickness(L-1))THEN
              IF(Root1stDepz_pft(N,NR,NZ).LE.CumSoilThickness(L))THEN
                RTDPP=Root1stDepz_pft(N,NR,NZ)+CanPHeight4WatUptake(NZ)
                Root1stSink_pvr(N,L,NR)=RTSK(iPlantRootProfile_pft(NZ))&
                  *RootAreaPopu*Root1stRadius_pvr(N,L,NZ)**2._r8/RTDPP
                RootSinkC(N)=RootSinkC(N)+Root1stSink_pvr(N,L,NR)
                RootSinkC_vr(N,L)=RootSinkC_vr(N,L)+Root1stSink_pvr(N,L,NR)
              ENDIF
            ENDIF

!
!     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
!     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
!
!     Root1stLocDepz_vr=depth of primary root axis in layer
!     RTDP1=primary root depth from soil surface
!     CumSoilThickness=depth from soil surface to layer bottom
!     RTDPX=distance behind growing point for secondary roots
!     DLYR=layer thickness
!     SeedDepth_pft=seeding depth
!     HypoctoHeight_pft=hypocotyledon height
!     CanPHeight4WatUptake=canopy height for water uptake
!     RTDPS=secondary root depth from canopy
!     RTSKP,RTSKS=primary,secondary root sink strength
!     RTN2=number of secondary root axes
!     Root2ndSink_pvr=total secondary root sink strength
!     AveLen2ndRoot_pvr=average secondary root length
!     RootSinkC,RootSinkC_vr=total root sink strength
!
            Root1stLocDepz_vr(NR,L)=AZMAX1(Root1stDepz_pft(ipltroot,NR,NZ)-CumSoilThickness(L-1)-RTDPX)
            Root1stLocDepz_vr(NR,L)=AZMAX1(AMIN1(DLYR3(L),Root1stLocDepz_vr(NR,L)) &
              -AZMAX1(SeedDepth_pft(NZ)-CumSoilThickness(L-1)-HypoctoHeight_pft(NZ)))
            RTDPS=AMAX1(SeedDepth_pft(NZ),CumSoilThickness(L-1))+0.5_r8*Root1stLocDepz_vr(NR,L) &
              +CanPHeight4WatUptake(NZ)
            IF(RTDPS.GT.ZERO)THEN
              RTSKP=RootAreaPopu*Root1stRadius_pvr(N,L,NZ)**2._r8/RTDPS
              RTSKS=safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)*Radius2ndRoot_pvr(N,L,NZ)**2._r8,AveLen2ndRoot_pvr(N,L,NZ))
              IF(RTSKP+RTSKS.GT.ZEROP(NZ))THEN
                Root2ndSink_pvr(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
              ELSE
                Root2ndSink_pvr(N,L,NR)=0._r8
              ENDIF
            ELSE
              Root2ndSink_pvr(N,L,NR)=0._r8
            ENDIF
          ELSE
            Root2ndSink_pvr(N,L,NR)=safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)&
              *Radius2ndRoot_pvr(N,L,NZ)**2._r8,AveLen2ndRoot_pvr(N,L,NZ))
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
