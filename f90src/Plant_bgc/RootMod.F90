module RootMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod  , only : safe_adb,AZMAX1,AZMIN1
  use EcoSIMCtrlMod, only : lverb  
  use EcosimConst
  use GrosubPars
  use ElmIDMod
  use PlantMathFuncMod
  use PlantAPIData
  use NoduleBGCMod
  use PlantBalMod, only : SumRootBiome
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: RootBGCModel
  contains

  subroutine RootBGCModel(I,J,NZ,BegRemoblize,PTRT,TFN6_vr,CNRTW,CPRTW,RootPrimeAxsNum,&
    RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)
!
! Roots and mycorrhizae both have primary and secondary components
! primary roots are accordant to shoots, seminar root corresponds to main branch.
! not all plants have seminal roots.
  implicit none
  integer , intent(in) :: I,J,NZ
  integer , intent(in) :: BegRemoblize  !remobilization flag
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: CNRTW,CPRTW
  real(r8), intent(in) :: RootPrimeAxsNum
  real(r8), intent(in) :: PTRT      !shoot-root nonstrucal C/N exchange modifier
  real(r8), intent(out) :: RootSinkC_vr(jroots,JZ1)
  real(r8), intent(out) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(out) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(out) :: RootSinkC(jroots)

  integer, parameter  :: NumRootAxes4DeadPlant =0    !
  real(r8) :: TotRootVol

  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)

  associate(                                                      &
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft,    &
    ZERO4LeafVar_pft        => plt_biom%ZERO4LeafVar_pft,         &
    DLYR3                   => plt_site%DLYR3,                    &
    PlantPopulation_pft     => plt_site%PlantPopulation_pft,      &
    ZERO                    => plt_site%ZERO,                     &
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft, &
    iPlantRootState_pft     => plt_pheno%iPlantRootState_pft,     &
    iPlantShootState_pft    => plt_pheno%iPlantShootState_pft,    &
    fRootGrowPSISense_pvr   => plt_pheno%fRootGrowPSISense_pvr,   &
    SeedMeanLen_pft         => plt_morph%SeedMeanLen_pft,         &
    RootVH2O_pvr            => plt_morph%RootVH2O_pvr,            &
    RootAreaPerPlant_pvr    => plt_morph%RootAreaPerPlant_pvr,    &
    RootPoreVol_pvr         => plt_morph%RootPoreVol_pvr,         &
    RootLenDensPerPlant_pvr => plt_morph%RootLenDensPerPlant_pvr, &
    RootLenPerPlant_pvr     => plt_morph%RootLenPerPlant_pvr,     &
    NGTopRootLayer_pft      => plt_morph%NGTopRootLayer_pft,      &
    RootPorosity_pft        => plt_morph%RootPorosity_pft,        &
    SeedVolumeMean_pft      => plt_morph%SeedVolumeMean_pft,      &
    SeedAreaMean_pft        => plt_morph%SeedAreaMean_pft,        &
    NumRootAxes_pft         => plt_morph%NumRootAxes_pft          &
  )
!     ROOT GROWTH
!
  
!  call RootCheck(I,J,NZ,'head')
  if(lverb)write(*,*)'SummarizeRootSink'
  call SummarizeRootSink(I,J,NZ,RootPrimeAxsNum,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)

  if(lverb)write(*,*)'RootBiochemistry'
  call RootBiochemistry(I,J,NZ,TFN6_vr,CNRTW,CPRTW,RootPrimeAxsNum,RootSinkC_vr,&
    Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)
 
!  call RootCheck(I,J,NZ,'rootbiochem')
    
!
!     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
!     GERMINATION)
!
!   call SumRootBiome(NZ,mass_inital)

  RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedMeanLen_pft(NZ)
  IF(DLYR3(NGTopRootLayer_pft(NZ)).GT.ZERO)THEN
    !by m3
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)= &
      RootLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)/DLYR3(NGTopRootLayer_pft(NZ))
  ELSE
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=0._r8
  ENDIF
  TotRootVol=RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedVolumeMean_pft(NZ)*PlantPopulation_pft(NZ)
  RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)      = RootPorosity_pft(ipltroot,NZ)*TotRootVol
  RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)         = (1.0_r8-RootPorosity_pft(ipltroot,NZ))*TotRootVol
  RootAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+&
    SeedAreaMean_pft(NZ)

  IF(NumRootAxes4DeadPlant.EQ.NumRootAxes_pft(NZ) .OR. (SeasonalNonstElms_pft(ielmc,NZ).LE.ZERO4LeafVar_pft(NZ).AND. &
    iPlantPhenolPattern_pft(NZ).NE.iplt_annual))THEN
    iPlantRootState_pft(NZ)=iDead
    iPlantShootState_pft(NZ)=iDead
  ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
  call RootNodulBiochemistry(I,J,NZ,TFN6_vr)

!   call SumRootBiome(NZ,mass_finale)
  
!  call RootCheck(I,J,NZ,'tail')
  end associate
  end subroutine RootBGCModel

!------------------------------------------------------------------------------------------

  subroutine RootBiochemistry(I,J,NZ,TFN6_vr,CNRTW,CPRTW,RootPrimeAxsNum,&
    RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1),CNRTW,CPRTW
  real(r8), intent(in) :: RootPrimeAxsNum
  REAL(R8), INTENT(in) :: RootSinkC_vr(jroots,JZ1)
  real(r8), INTENT(in) :: RootSinkC(jroots)
  real(r8), INTENT(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)

  integer :: LL,LZ,L1,L,K,lx,M,NR,N,NTG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: NonstElmGradt
  real(r8) :: CPOOLX
  real(r8) :: FWTRT
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
  integer :: iRootXsUpdateFlag(jroots,JZ1)  
  INTEGER :: NRX(jroots,JZ1)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: masst_inital(NumPlantChemElms)
  real(r8) :: masst_finale(NumPlantChemElms)  
  real(r8) :: litrflx(NumPlantChemElms)
  real(r8) :: RCO2flx
  real(r8) :: litrflxt(NumPlantChemElms)
  real(r8) :: RCO2flxt
  character(len=120) :: procmsg

!     begin_execution
  associate(                                                             &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,        &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr,       &
    RootElms_pft              => plt_biom%RootElms_pft,                  &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,                &
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft,         &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,        &
    RootBiomGrosYld_pft       => plt_allom%RootBiomGrosYld_pft,          &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,        &
    fRootGrowPSISense_pvr     => plt_pheno%fRootGrowPSISense_pvr ,       &
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                     &
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,                  &
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft,        &
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr,               &
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr,               &
    SoilResit4RootPentrate_vr => plt_soilchem%SoilResit4RootPentrate_vr, &
    VLSoilPoreMicP_vr         => plt_soilchem%VLSoilPoreMicP_vr,         &
    NU                        => plt_site%NU,                            &
    ZERO                      => plt_site%ZERO,                          &
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,           &
    ZEROS2                    => plt_site%ZEROS2,                        &
    DLYR3                     => plt_site%DLYR3,                         &
    NL                        => plt_site%NL,                            &
    k_fine_litr               => pltpar%k_fine_litr,                     &
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr,              &
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr,      &
    RootAreaPerPlant_pvr      => plt_morph%RootAreaPerPlant_pvr,         &
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr,                 &
    Root1stMaxRadius1_pft     => plt_morph%Root1stMaxRadius1_pft,        &
    Root2ndMaxRadius1_pft     => plt_morph%Root2ndMaxRadius1_pft,        &
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft,         &
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr,            &
    RootLenPerPlant_pvr       => plt_morph%RootLenPerPlant_pvr,          &
    Root2ndAveLen_pvr         => plt_morph%Root2ndAveLen_pvr,            &
    Root2ndRadius_pvr         => plt_morph%Root2ndRadius_pvr,            &
    Root2ndXSecArea_pft       => plt_morph%Root2ndXSecArea_pft,          &
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft,         &
    Root2ndXNum_pvr           => plt_morph%Root2ndXNum_pvr,              &
    Root1stXSecArea_pft       => plt_morph%Root1stXSecArea_pft,          &
    RootPorosity_pft          => plt_morph%RootPorosity_pft,             &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,             &
    RootVolPerMassC_pft       => plt_morph%RootVolPerMassC_pft,          &
    SeedMeanLen_pft           => plt_morph%SeedMeanLen_pft,              &
    MY                        => plt_morph%MY,                           &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,           &
    NIXBotRootLayer_pft       => plt_morph%NIXBotRootLayer_pft           &
  )

  iRootXsUpdateFlag       = ifalse
  NRX                     = ifalse
  NIXBotRootLayer_pft(NZ) = NGTopRootLayer_pft(NZ)
!
!     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
!
!  mass_inital(ielmc)=mass_inital(ielmc)+SeasonalNonstElms_pft(ielmc,NZ)
  litrflx=0._r8;RCO2flx=0._r8
  D5010: DO N=1,MY(NZ)
    D5000: DO L=NU,MaxSoiL4Root_pft(NZ)
!
!     IDENTIFY NEXT LOWER ROOT LAYER
!
!     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
!
      IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
        !why not using MaxSoiL4Root_pft(NZ)
        D5003: DO LZ=L+1,NL
          IF(VLSoilPoreMicP_vr(LZ).GT.ZEROS2 .OR. LZ.EQ.NL)THEN
            L1=LZ
            EXIT
          ENDIF
        ENDDO D5003

        litrflxt=0._r8;RCO2flxt=0._r8
      !   call SumRootBiome(NZ,masst_inital)

        IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
          fRootGrowPSISense_pvr(N,L,NZ)=EXP(0.05_r8*AMAX1(PSIRoot_pvr(N,L,NZ),-5000._r8))
        ELSE
          fRootGrowPSISense_pvr(N,L,NZ)=EXP(0.10_r8*AMAX1(PSIRoot_pvr(N,L,NZ),-5000._r8))
        ENDIF  
        
        call GrowRootMycoAxes(I,J,N,L,L1,NZ,NRX,iRootXsUpdateFlag,TFN6_vr,&
          RootPrimeAxsNum,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,CNRTW,CPRTW,&
          fRootGrowPSISense_pvr(N,L,NZ),TotRoot2ndLen,TotRoot1stLen,Root2ndC,Root1stC,litrflxt,RCO2flxt)
        
      !   call SumRootBiome(NZ,masst_finale)
        litrflx=litrflx+litrflxt
        RCO2flx=RCO2flx+RCO2flxt

!====================================================================================
        masst_inital=masst_finale
!     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
!     SEASONAL STORAGE POOL IS DEPLETED
!
!     RootMycoActiveBiomC_pvr,WTRT=total root C mass
!     WTRVC=storage C
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!     CPOOLR=non-structural C mass in root
!     Note, 03/28/2024, Jinyun Tang: I need to take a careful look at the following code, particularly for WTRTLX and WTRTTX,
!     which are essentially the same. However, assuming the flux is driven by C concentration of the seasonal storage
!     and the nonstructural biomass in root, then WTRTTX should be defined with stalk volume.
!     Question: where is seasonal storage located?
        IF(L.LE.NIXBotRootLayer_pft(NZ))THEN
          IF(RootMycoActiveBiomC_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ) .AND. RootElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) &
            .AND. SeasonalNonstElms_pft(ielmc,NZ).LT.XFRX*RootElms_pft(ielmc,NZ))THEN
            FWTRT                                = RootMycoActiveBiomC_pvr(N,L,NZ)/RootElms_pft(ielmc,NZ)
            WTRTLX                               = RootMycoActiveBiomC_pvr(N,L,NZ)
            WTRTTX                               = RootElms_pft(ielmc,NZ)*FWTRT
            WTRTTT                               = WTRTLX+WTRTTX
            CPOOLX                               = AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
            WTRVCX                               = AZMAX1(SeasonalNonstElms_pft(ielmc,NZ)*FWTRT)
            NonstElmGradt                        = (WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
            XFRC                                 = XFRY*AZMIN1(NonstElmGradt)
            RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+XFRC
            SeasonalNonstElms_pft(ielmc,NZ)      = SeasonalNonstElms_pft(ielmc,NZ)-XFRC
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
!     RRAD1,Root2ndRadius_pvr=primary,secondary root radius
!     RootAreaPerPlant_pvr=root surface area per plant
!     Root2ndAveLen_pvr=average secondary root length
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!
        IF(N.EQ.ipltroot)THEN
          TotRoot1stLen=TotRoot1stLen*FracRootElmAlloc2Litr(ielmc,k_fine_litr)
          TotRoot2ndLen=TotRoot2ndLen*FracRootElmAlloc2Litr(ielmc,k_fine_litr)
        ENDIF

        TotPopuRoot1stLen_rpvr = TotRoot1stLen*PlantPopulation_pft(NZ)
        TotPopuRootLen         = TotRoot2ndLen+TotPopuRoot1stLen_rpvr
        WTRTT                  = Root2ndC+Root1stC
        IF(TotPopuRootLen.GT.ZERO4Groth_pft(NZ) .AND. WTRTT.GT.ZERO4Groth_pft(NZ) &
          .AND. PlantPopulation_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
          RootLenPerPlant_pvr(N,L,NZ)=TotPopuRootLen/PlantPopulation_pft(NZ)
          IF(DLYR3(L).GT.ZERO)THEN
            !per volume
            RootLenDensPerPlant_pvr(N,L,NZ)=RootLenPerPlant_pvr(N,L,NZ)/DLYR3(L)
          ELSE
            RootLenDensPerPlant_pvr(N,L,NZ)=0._r8
          ENDIF
          TotRootVol=AMAX1(Root1stXSecArea_pft(N,NZ)*TotPopuRoot1stLen_rpvr+Root2ndXSecArea_pft(N,NZ)*TotRoot2ndLen &
            ,WTRTT*RootVolPerMassC_pft(N,NZ)*PSIRootTurg_vr(N,L,NZ))
          RootPoreVol_pvr(N,L,NZ) = RootPorosity_pft(N,NZ)*TotRootVol
          RootVH2O_pvr(N,L,NZ)    = (1.0_r8-RootPorosity_pft(N,NZ))*TotRootVol
          !primary roots
          Root1stRadius_pvr(N,L,NZ)=AMAX1(Root1stMaxRadius1_pft(N,NZ),&
            (1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
          !secondary roots
          Root2ndRadius_pvr(N,L,NZ)=AMAX1(Root2ndMaxRadius1_pft(N,NZ),&
            (1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root2ndMaxRadius_pft(N,NZ))
          TotRootArea=TwoPiCON*(Root1stRadius_pvr(N,L,NZ)*TotPopuRoot1stLen_rpvr+Root2ndRadius_pvr(N,L,NZ)*TotRoot2ndLen)
          IF(Root2ndXNum_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
            Root2ndAveLen_pvr(N,L,NZ)=AMAX1(Root2ndAveLenMin,TotRoot2ndLen/Root2ndXNum_pvr(N,L,NZ))
          ELSE
            Root2ndAveLen_pvr(N,L,NZ)=Root2ndAveLenMin
          ENDIF
          RootAreaPerPlant_pvr(N,L,NZ)=TotRootArea/PlantPopulation_pft(NZ)
        ELSE
          RootLenPerPlant_pvr(N,L,NZ)     = 0._r8
          RootLenDensPerPlant_pvr(N,L,NZ) = 0._r8
          RootPoreVol_pvr(N,L,NZ)         = 0._r8
          RootVH2O_pvr(N,L,NZ)            = 0._r8
          RootAreaPerPlant_pvr(N,L,NZ)    = 0._r8
          Root1stRadius_pvr(N,L,NZ)       = Root1stMaxRadius_pft(N,NZ)
          Root2ndRadius_pvr(N,L,NZ)       = Root2ndMaxRadius_pft(N,NZ)
          Root2ndAveLen_pvr(N,L,NZ)       = Root2ndAveLenMin
          DO NTG=idg_beg,idg_NH3
            RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ) &
              -(trcg_rootml_pvr(NTG,N,L,NZ)+trcs_rootml_pvr(NTG,N,L,NZ))
          ENDDO
          trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
          trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
        ENDIF
      !   call SumRootBiome(NZ,masst_finale)
      ENDIF
    ENDDO D5000
  ENDDO D5010

!   call SumRootBiome(NZ,mass_finale)
!  mass_finale(ielmc)=mass_finale(ielmc)+SeasonalNonstElms_pft(ielmc,NZ)

  end associate
  end subroutine RootBiochemistry
!------------------------------------------------------------------------------------------

  subroutine Grow2ndRootAxes(I,J,L,NZ,N,NR,CNRTW,CPRTW,WFNR,WFNRG,DMRTD,fRootGrowPSISense,RootPrimeAxsNum,&
    TFN6_vr,RootSinkC_vr,Root2ndSink_pvr,litrflx2,RCO2flx2,Root2ndStrutRemob)
  implicit none
  integer, intent(in) :: I,J,L,NZ,N,NR
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: WFNR,WFNRG
  real(r8), intent(in) :: DMRTD    
  real(r8), intent(in) :: fRootGrowPSISense  
  real(r8), intent(in) :: RootPrimeAxsNum
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  REAL(R8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)  
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(out) :: litrflx2(NumPlantChemElms),RCO2flx2
  real(r8), intent(out) :: Root2ndStrutRemob(NumPlantChemElms)    
  real(r8) :: RCCC,RCCN,RCCP,FracRemobl  
  real(r8) :: Rmaint2nd_CO2,Frac2Senes2    
  real(r8) :: RNonstCO2_OUltd,RNonstCO2_Oltd  
  REAL(R8) :: RCO2Nonst4Xmaint2nd_Oltd,RCO2Nonst4Xmaint2nd_OUltd
  real(r8) :: RCO2Nonst4Nassim_Oltd,RCO2Nonst4Nassim_OUltd
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4GrowC_OUltd
  real(r8) :: RootMycoNonst4GrowC_Oltd    
  real(r8) :: Root2ndExtension  
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms)
  real(r8) :: Root2ndNetGrowthElms(NumPlantChemElms)
  real(r8) :: Remobl2ndcycl(NumPlantChemElms)
  real(r8) :: Remobl2ndelm(NumPlantChemElms)
  real(r8) :: dmass(NumPlantChemElms)  
  real(r8) :: RCO2XMaint2nd_Oltd
  real(r8) :: RGrowCO2_Oltd,RCO2T2nd_Oltd
  real(r8) :: RCO2XMaint2nd_OUltd
  real(r8) :: RGrowCO2_OUltd
  real(r8) :: RCO2T2nd_OUltd
  real(r8) :: CNPG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: DMRTR
  real(r8) :: PPOOLB,ZPOOLB
  real(r8) :: RTN2X,RTN2Y  
  real(r8) :: dsenecE,dfineE,dwoodyE
  real(r8) :: respscal  
  real(r8) :: FNP
  real(r8) :: FRTN  
  integer  :: M,NE

  associate(                                                             &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr,     &
    RootNonstructElmConc_rpvr  => plt_biom%RootNonstructElmConc_rpvr,      &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,        &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,              &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr,       &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,                &
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                     &
    RootCO2Autor_pvr          => plt_rbgc%RootCO2Autor_pvr,              &
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr,            &
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr,            &
    RAutoRootO2Limter_rpvr     => plt_rbgc%RAutoRootO2Limter_rpvr,         &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,          &
    rCNNonstRemob_pft         => plt_allom%rCNNonstRemob_pft,            &
    rCPNonstRemob_pft         => plt_allom%rCPNonstRemob_pft,            &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,        &
    CNRTS_pft                 => plt_allom%CNRTS_pft,                    &
    CPRTS_pft                 => plt_allom%CPRTS_pft,                    &
    RootBiomGrosYld_pft       => plt_allom%RootBiomGrosYld_pft,          &
    k_woody_litr              => pltpar%k_woody_litr,                    &
    k_fine_litr               => pltpar%k_fine_litr,                     &
    icwood                    => pltpar%icwood,                          &
    iroot                     => pltpar%iroot,                           &
    inonstruct                => pltpar%inonstruct,                      &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,        &
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft,         &
    fTgrowRootP_vr            => plt_pheno%fTgrowRootP_vr,               &
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch,          &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,          &
    DLYR3                     => plt_site%DLYR3,                         &
    ZERO                      => plt_site%ZERO,                          &
    Root2ndXNum_pvr           => plt_morph%Root2ndXNum_pvr,              &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,           &
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr,             &
    RootBranchFreq_pft        => plt_morph%RootBranchFreq_pft,           &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,                &
    Root2ndSpecLen_pft        => plt_morph%Root2ndSpecLen_pft,           &
    Root2ndLen_rpvr            => plt_morph%Root2ndLen_rpvr,               &
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft,            &
    C4PhotosynDowreg_brch     => plt_photo%C4PhotosynDowreg_brch         &
  )
!     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
!     TO CURRENT AXIS
!
!     Root2ndSink_pvr=total secondary root sink strength
!     RootSinkC_vr=total root sink strength
!     FRTN=fraction of secondary root sink strength in axis
!
  dmass=0._r8;litrflx2=0._r8
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
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
  IF(RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CNPG=AMIN1(RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmn,N,L,NZ) &
      +RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CNKI),RootNonstructElmConc_rpvr(ielmp,N,L,NZ) &
      /(RootNonstructElmConc_rpvr(ielmp,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CPKI))
  ELSE
    CNPG=1.0_r8
  ENDIF
!
!     SECONDARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     Rmaint2nd_CO2=root maintenance respiration
!     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT2N=secondary root N mass
!     TFN6_vr=temperature function for root maintenance respiration
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     fRootGrowPSISense=growth function of root water potential
!
  Rmaint2nd_CO2=AZMAX1(RmSpecPlant*RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ))*TFN6_vr(L)
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)) .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint2nd_CO2=Rmaint2nd_CO2*fRootGrowPSISense
  ENDIF
!
!     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RNonstCO2_OUltd=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     
!     FRTN=fraction of secondary root sink strength in axis
!     CPOOL=non-structural C mass
!     fTgrowRootP_vr=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!     fRootGrowPSISense=growth function of root water potential
!
  respscal=VMXC*FRTN*fTgrowRootP_vr(L,NZ)*fRootGrowPSISense &
    *CNPG*C4PhotosynDowreg_brch(MainBranchNum_pft(NZ),NZ)

  RNonstCO2_OUltd=AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))*respscal

!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!     RCO2XMaint2nd_OUltd,RCO2XMaint2nd_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
  RNonstCO2_Oltd      = RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RCO2XMaint2nd_OUltd = RNonstCO2_OUltd-Rmaint2nd_CO2
  RCO2XMaint2nd_Oltd  = RNonstCO2_Oltd-Rmaint2nd_CO2
  RGrowCO2_OUltd      = AZMAX1(RCO2XMaint2nd_OUltd)*WFNRG
  RGrowCO2_Oltd       = AZMAX1(RCO2XMaint2nd_Oltd)*WFNRG
!
!     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS_pft,CPRTS_pft=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!
  DMRTR  = DMRTD*FRTN
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))
  FNP    = AMIN1(ZPOOLB/CNRTS_pft(NZ),PPOOLB/CPRTS_pft(NZ))*DMRTR
  IF(RGrowCO2_OUltd.GT.0.0_r8)THEN
    RGrowCO2_OUltd=AMIN1(RGrowCO2_OUltd,FNP)
  ELSE
    RGrowCO2_OUltd=0._r8
  ENDIF
  !active growth
  IF(RGrowCO2_Oltd.GT.0.0)THEN
    RGrowCO2_Oltd=AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ))
  ELSE
    RGrowCO2_Oltd=0._r8
  ENDIF

!
!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN 'READQ'
!
!     RootMycoNonst4GrowC_OUltd,RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn unltd,ltd by O2
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption, i.e. respiration quotient
!     RootMycoNonst4Grow_OUltd(ielmc),RootMycoNonst4Grow_Oltd(ielmc)=root C growth unltd,ltd by O2
!     RootBiomGrosYld_pft=root growth yield
!     RootMycoNonst4Grow_OUltd(:)=nonstructural N,P unlimited,limited by O2 used in growth
!     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
!
  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRTD
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRTD
  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)
  RCO2Nonst4Nassim_OUltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_OUltd(ielmn))

  RootMycoNonst4Grow_Oltd(ielmc)=RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)
!  print*,'rootmyco4gro',RootMycoNonst4Grow_Oltd(ielmc),RootMycoNonst4GrowC_Oltd,RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmn,N,L,NZ) &
    ,RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
  RootMycoNonst4Grow_Oltd(ielmp)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmp,N,L,NZ) &
    ,RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))
  RCO2Nonst4Nassim_Oltd=AZMAX1(1.70_r8*RootMycoNonst4Grow_Oltd(ielmn))
!
!     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LitrFall
!
!     iPlantCalendar_brch(ipltcal_Emerge),=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0.AND.&
    RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmn,N,L,NZ)&
      ,RootNonstructElmConc_rpvr(ielmn,N,L,NZ) &
      +RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CNKI) &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmp,N,L,NZ),RootNonstructElmConc_rpvr(ielmp,N,L,NZ) &
      +RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8 &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ),RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+ &
      RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8 &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ),RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+ &
      RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/CPKI)))
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
!     RCO2XMaint2nd_OUltd,RCO2XMaint2nd_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2Nonst4Xmaint2nd_OUltd,RCO2Nonst4Xmaint2nd_Oltd=excess maintenance respiration unltd,ltd by O2
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!     Root2ndStrutRemob(:)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     Frac2Senes2=fraction of secondary root C to be remobilized
!
  IF(-RCO2XMaint2nd_OUltd.GT.0.0_r8)THEN
    IF(-RCO2XMaint2nd_OUltd.LT.RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
      RCO2Nonst4Xmaint2nd_OUltd=-RCO2XMaint2nd_OUltd
    ELSE
      RCO2Nonst4Xmaint2nd_OUltd=AZMAX1(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)
    ENDIF
  ELSE
    RCO2Nonst4Xmaint2nd_OUltd=0._r8
  ENDIF

  !there is maintenance deficit,
  IF(-RCO2XMaint2nd_Oltd.GT.0.0_r8)THEN
    IF(-RCO2XMaint2nd_Oltd.LT.RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
      RCO2Nonst4Xmaint2nd_Oltd=-RCO2XMaint2nd_Oltd
    ELSE
      RCO2Nonst4Xmaint2nd_Oltd=AZMAX1(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)*RAutoRootO2Limter_rpvr(N,L,NZ)
    ENDIF
  ELSE
    RCO2Nonst4Xmaint2nd_Oltd=0._r8
  ENDIF

  !Maintenance deficit leads to remobilization/degradation of structural biomass to pay off the deficit
  IF(RCO2Nonst4Xmaint2nd_Oltd.GT.0.0_r8 .AND. RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    Root2ndStrutRemob(ielmc)=RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC
    FracRemobl=Root2ndStrutRemob(ielmc)/RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
    Root2ndStrutRemob(ielmn)=RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ)*(RCCN+(1.0_r8-RCCN)*FracRemobl)
    Root2ndStrutRemob(ielmp)=RootMyco2ndStrutElms_rpvr(ielmp,N,L,NR,NZ)*(RCCP+(1.0_r8-RCCP)*FracRemobl)

    IF(Root2ndStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
      Frac2Senes2=AZMAX1(AMIN1(1.0_r8,RCO2Nonst4Xmaint2nd_Oltd/Root2ndStrutRemob(ielmc)))
    ELSE
      Frac2Senes2=1.0_r8
    ENDIF
  ELSE
    Root2ndStrutRemob(1:NumPlantChemElms)=0._r8
    Frac2Senes2=0._r8
  ENDIF
!
!     SECONDARY ROOT LitrFall CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     Frac2Senes2=fraction of secondary root C to be remobilized
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     Root2ndStrutRemob(ielmc),Root2ndStrutRemob(ielmn),Root2ndStrutRemob(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
  if(Frac2Senes2>0._r8)then 
    DO NE=1,NumPlantChemElms
      Remobl2ndelm(NE)  = RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*Frac2Senes2
      Remobl2ndcycl(NE) = Root2ndStrutRemob(NE)*Frac2Senes2
    ENDDO  
    D6350: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        dsenecE                                      = AZMAX1(Remobl2ndelm(NE)-Remobl2ndcycl(NE))
        dwoodyE                                      = ElmAllocmat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_woody_litr)
        dfineE                                       = ElmAllocmat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_fine_litr)
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ) = LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)+dwoodyE
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)  = LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+dfineE
        litrflx2(NE)                                 = litrflx2(NE)+dwoodyE+dfineE
        dmass(NE)                                    = dmass(NE)-dfineE-dwoodyE
      ENDDO    
    ENDDO D6350
  endif
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     Rmaint2nd_CO2=root maintenance respiration
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn ltd by O2
!     RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
!     RCO2Nonst4Xmaint2nd_Oltd=excess maintenance respiration ltd by O2
!     Frac2Senes2=fraction of secondary root C to be remobilized
!     Root2ndStrutRemob(:)=remobilization of C,N,P from senescing root
!     RootMycoNonst4Grow_Oltd(:)=nonstructural N,P ltd by O2 used in growth
!
    DO NE=1,NumPlantChemElms
!          print*,'conssss',RootMycoNonstElms_rpvr(NE,N,L,NZ),RootMycoNonst4Grow_Oltd(NE)
      RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Grow_Oltd(NE)
      dmass(NE)=dmass(NE)+RootMycoNonst4Grow_Oltd(NE)

    ENDDO

    if(Frac2Senes2>0._r8)then
!      Remobilization is added to nonstrucal metabolites       
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+Remobl2ndcycl(NE)
      dmass(NE)=dmass(NE)-Remobl2ndcycl(NE)
    ENDDO
    endif

!
!     TOTAL SECONDARY ROOT RESPIRATION
!
!     RCO2T2nd_OUltd,RCO2T2nd_Oltd=total C respiration unltd,ltd by O2
!     Rmaint2nd_CO2=root maintenance respiration
!     RNonstCO2_OUltd,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     RCO2Nonst4Xmaint2nd_OUltd,RCO2Nonst4Xmaint2nd_Oltd=excess maintenance respiration unltd,ltd by O2
!     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
!     RootCO2Autor_pvr=total root respiration
!     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
!
  RCO2T2nd_OUltd                       = AMIN1(Rmaint2nd_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd+RCO2Nonst4Nassim_OUltd+RCO2Nonst4Xmaint2nd_OUltd
  RCO2T2nd_Oltd                        = AMIN1(Rmaint2nd_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+RCO2Nonst4Nassim_Oltd+RCO2Nonst4Xmaint2nd_Oltd
  RCO2T2nd_Oltd                        = AMIN1(RCO2T2nd_Oltd,AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)))
  RootRespPotent_pvr(N,L,NZ)           = RootRespPotent_pvr(N,L,NZ)+RCO2T2nd_OUltd
  RootCO2EmisPot_pvr(N,L,NZ)           = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T2nd_Oltd
  RootCO2Autor_pvr(N,L,NZ)             = RootCO2Autor_pvr(N,L,NZ)-RCO2T2nd_Oltd
  RCO2flx2                             = -RCO2T2nd_Oltd
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T2nd_Oltd

!
!     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     Root2ndExtension=secondary root length extension
!     RootMycoNonst4Grow_Oltd(ielmc)=secondary root C growth ltd by O2
!     Root2ndSpecLen_pft=specific secondary root length from startq.f
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     Frac2Senes2=fraction of secondary root C to be remobilized
!     Root2ndLen_rpvr=secondary root length
!     Root2ndNetGrowthElms(ielmc),Root2ndNetGrowthElms(ielmn),Root2ndNetGrowthElms(ielmp)=net root C,N,P growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RootMycoNonst4Grow_Oltd(ielmn),RootMycoNonst4Grow_Oltd(ielmp)=nonstructural N,P ltd by O2 used in growth
!
  Root2ndExtension=RootMycoNonst4Grow_Oltd(ielmc)*Root2ndSpecLen_pft(N,NZ)*WFNR*FracRootElmAlloc2Litr(ielmc,k_fine_litr) 
  Root2ndNetGrowthElms=RootMycoNonst4Grow_Oltd

  if(Frac2Senes2>0._r8)then
    Root2ndExtension=Root2ndExtension-Frac2Senes2*Root2ndLen_rpvr(N,L,NR,NZ)
    DO NE=1,NumPlantChemElms
      Root2ndNetGrowthElms(NE) = Root2ndNetGrowthElms(NE)-Remobl2ndelm(NE)
      dmass(NE)                = dmass(NE)+Remobl2ndelm(NE)
    ENDDO
  endif
!
!     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
!     AND AXIS NUMBER
!
!     Root2ndLen_rpvr=secondary root length
!     Root2ndExtension=secondary root length extension
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     Root2ndNetGrowthElms(ielmc),Root2ndNetGrowthElms(ielmn),Root2ndNetGrowthElms(ielmp)=net root C,N,P growth
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
!     RootBranchFreq_pft=root branching frequency from PFT file
!     Root2ndXNum_rpvr,Root2ndXNum_pvr=number of secondary root axes
!
  Root2ndLen_rpvr(N,L,NR,NZ)=Root2ndLen_rpvr(N,L,NR,NZ)+Root2ndExtension
  DO NE=1,NumPlantChemElms
    RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)+Root2ndNetGrowthElms(NE)
  ENDDO
  RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)+AMIN1(rCNNonstRemob_pft(NZ) &
    *RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ) &
    ,rCPNonstRemob_pft(NZ)*RootMyco2ndStrutElms_rpvr(ielmp,N,L,NR,NZ))

  !secondary root axes addition is a quadratic function of branching frequency
  RTN2X                       = RootBranchFreq_pft(NZ)*RootPrimeAxsNum
  RTN2Y                       = RootBranchFreq_pft(NZ)*RTN2X
  Root2ndXNum_rpvr(N,L,NR,NZ) = (RTN2X+RTN2Y)*DLYR3(L)
  Root2ndXNum_pvr(N,L,NZ)     = Root2ndXNum_pvr(N,L,NZ)+Root2ndXNum_rpvr(N,L,NR,NZ)
  end associate  
  end subroutine Grow2ndRootAxes  

!------------------------------------------------------------------------------------------

  subroutine GrowRootMycoAxes(I,J,N,L,L1,NZ,NRX,iRootXsUpdateFlag,TFN6_vr,&
    RootPrimeAxsNum,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,CNRTW,CPRTW,&
    fRootGrowPSISense,TotRoot2ndLen,TotRoot1stLen,Root2ndC,Root1stC,litrflx,RCO2flx)
  !
  !Description
  !Grow root axes in laye L  
  implicit none
  integer , intent(in) :: I,J
  INTEGER , INTENT(IN) :: N      !root type id
  integer , intent(in) :: L      !current root layer
  integer , intent(in) :: L1     !next root growable layer, soil porosity >0 
  integer , intent(in) :: NZ     !pft id
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: RootPrimeAxsNum
  REAL(R8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: CNRTW,CPRTW
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(out) :: TotRoot2ndLen,TotRoot1stLen,Root2ndC,Root1stC
  integer , intent(inout) :: NRX(jroots,JZ1)  
  integer , intent(inout) :: iRootXsUpdateFlag(jroots,JZ1)  
  real(r8), intent(out) :: litrflx(NumPlantChemElms)
  real(r8), intent(out) :: RCO2flx
  real(r8) :: Root2ndStrutRemob(NumPlantChemElms)    
  real(r8) :: DMRTD  
  real(r8) :: WFNR,WFNRG  
  real(r8) :: SoilRest4Root2ndPentrate  
  real(r8) :: FRCO2
  real(r8) :: FSNCM
  real(r8) :: FSNCP
  integer  :: NR
  real(r8) :: litrflxt(NumPlantChemElms),litrflx2(NumPlantChemElms)
  real(r8) :: RCO2flxt,RCO2flx2
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: massr1st(NumPlantChemElms),massr1st1(NumPlantChemElms)
  real(r8) :: massr2nd(NumPlantChemElms),massr2nd1(NumPlantChemElms)
  real(r8) :: massnonst(NumPlantChemElms),massnonst1(NumPlantChemElms)
  real(r8) :: massnodul(NumPlantChemElms),massnodul1(NumPlantChemElms)

!begin_execution
  associate(                                                             &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr,     &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,              &
    Root2ndRadius_pvr         => plt_morph%Root2ndRadius_pvr,            &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,        &
    NIXBotRootLayer_rpft      => plt_morph%NIXBotRootLayer_rpft,         &
    RootBiomGrosYld_pft       => plt_allom%RootBiomGrosYld_pft,          &
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr,              &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr,     &
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,                  &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,              &
    SoilResit4RootPentrate_vr => plt_soilchem%SoilResit4RootPentrate_vr, &
    NIXBotRootLayer_pft       => plt_morph%NIXBotRootLayer_pft           &
  )
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     SoilResit4RootPentrate_vr,SoilRest4Root2ndPentrate=soil resistance to secondary root penetration (MPa)
!     Root2ndRadius_pvr=secondary root radius
!     WFNR=water function for root extension
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     fRootGrowPSISense,WFNRG=growth,respiration function of root water potential
!     PSIRoot_pvr,PSIRootTurg_vr=root total,turgor water potential
!     DMRT=root growth yield
!
!   call SumRootBiome(NZ,mass_inital,massr1st1,massr2nd1,massnonst1,massnodul1)

  SoilRest4Root2ndPentrate = SoilResit4RootPentrate_vr(L)*Root2ndRadius_pvr(N,L,NZ)/1.0E-03_r8
  WFNR                     = AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-PSIMin4OrganExtens-SoilRest4Root2ndPentrate))
  WFNRG                    = fRespWatSens(WFNR,iPlantRootProfile_pft(NZ))


  !respiration fraction
  DMRTD         = 1.0_r8-RootBiomGrosYld_pft(NZ)
  TotRoot1stLen = 0._r8;TotRoot2ndLen = 0._r8
  Root1stC      = 0._r8;Root2ndC      = 0._r8
  litrflx       = 0._r8;RCO2flx       = 0._r8

!  print*
!
!     FOR EACH ROOT AXIS

  !=================================================================================
  !loop through all axes 
  D5050: DO NR=1,NumRootAxes_pft(NZ)
      
!     SECONDARY ROOT EXTENSION
!   make sure current layer is no deeper than the lowest root layer
    IF(L.LE.NIXBotRootLayer_rpft(NR,NZ) .AND. NRX(N,NR).EQ.ifalse)THEN
!
     call Grow2ndRootAxes(I,J,L,NZ,N,NR,CNRTW,CPRTW,WFNR,WFNRG,DMRTD,fRootGrowPSISense,RootPrimeAxsNum,&
       TFN6_vr,RootSinkC_vr,Root2ndSink_pvr,litrflx2,RCO2flx2,Root2ndStrutRemob)
     TotRoot2ndLen=TotRoot2ndLen+Root2ndLen_rpvr(N,L,NR,NZ)
     Root2ndC=Root2ndC+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
!     print*,'primroot',NR,L
!
!     PRIMARY ROOT EXTENSION
!
!     SoilBulkDensity_vr=soil bulk density
!     RTDP1,Root1stDepz2Surf=primary root depth from soil surface
!     CumSoilThickness_vr=depth from soil surface to layer bottom
!     ICHKL=flag for identifying layer with primary root tip
!     RTN1=number of primary root axes
!     RootPrimeAxsNum=multiplier for number of primary root axes
!
      IF(N.EQ.ipltroot)THEN

        call Grow1stRootAxes(I,J,N,NR,L,L1,NZ,RootPrimeAxsNum,Root1stSink_pvr,&
          Root2ndSink_pvr,RootSinkC_vr,Root2ndStrutRemob,fRootGrowPSISense,TFN6_vr,&
          DMRTD,CNRTW,CPRTW,TotRoot1stLen,Root1stC,iRootXsUpdateFlag,&
          NRX,litrflxt,RCO2flxt)
          litrflx=litrflx+litrflxt+litrflx2
          RCO2flx=RCO2flx+RCO2flxt+RCO2flx2

      ENDIF
!     TotRoot1stLen=total primary root length
!     Root1stC=total primary root C mass      
      TotRoot1stLen = TotRoot1stLen+Root1stLen_rpvr(N,L,NR,NZ)
      Root1stC      = Root1stC+RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
    ENDIF
    NIXBotRootLayer_pft(NZ)=MAX(NIXBotRootLayer_pft(NZ),NIXBotRootLayer_rpft(NR,NZ))
  !   call SumRootBiome(NZ,mass_finale,massr1st,massr2nd,massnonst,massnodul)

  ENDDO D5050
  end associate
  end subroutine GrowRootMycoAxes

!------------------------------------------------------------------------------------------

  subroutine Grow1stRootAxes(I,J,N,NR,L,L1,NZ,RootPrimeAxsNum,Root1stSink_pvr,&
    Root2ndSink_pvr,RootSinkC_vr,Root2ndStrutRemob,fRootGrowPSISense,TFN6_vr,&
    DMRTD,CNRTW,CPRTW,TotRoot1stLen,Root1stC,iRootXsUpdateFlag,&
    NRX,litrflxt,RCO2flxt)
  !
  !grow root axes
  implicit none
  integer,  intent(in) :: I,J,N,NR,L,L1,NZ
  real(r8), intent(in) :: RootPrimeAxsNum
  REAL(R8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: Root2ndStrutRemob(NumPlantChemElms)  
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRTD  
  integer,  intent(inout) :: iRootXsUpdateFlag(jroots,JZ1)
  integer,  intent(inout) :: NRX(jroots,JZ1)  
  real(r8), intent(inout) :: TotRoot1stLen
  real(r8), intent(inout) :: Root1stC
  real(r8), intent(out) :: litrflxt(NumPlantChemElms)
  real(r8), intent(out) :: RCO2flxt
  real(r8) :: Root1stDepz2Surf
  real(r8) :: FRTN,Frac2Senes1
  real(r8) :: GRTWTM  
  real(r8) :: WFNR,WFNRG
  real(r8) :: CNPG,TFRCO2,FRCO2  
  real(r8) :: Rmaint1st_CO2
  real(r8) :: RNonstCO2_OUltd,RGrowCO2_Oltd,RGrowCO2_OUltd
  real(r8) :: RNonstCO2_Oltd,RCO2XMaint1st_OUltd,RCO2XMaint1st_Oltd
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4GrowC_OUltd
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms),RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd
  real(r8) :: FSNCM,FSNCP  
  real(r8) :: dsenecE
  real(r8) :: RootMycoNonst4GrowC_Oltd
  real(r8) :: dRootMycoElms(NumPlantChemElms)  
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  integer  :: NE,M,LL
  
  associate(                                                         &
    ZERO                      => plt_site%ZERO,                      &
    RootCO2Autor_pvr          => plt_rbgc%RootCO2Autor_pvr,          &
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr,        &
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr,        &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft,     &
    fTgrowRootP_vr            => plt_pheno%fTgrowRootP_vr,           &
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr,    &
    MaxNumRootLays            => plt_site%MaxNumRootLays,            &
    CumSoilThickness_vr       => plt_site%CumSoilThickness_vr,       &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,       &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr,   &
    RootNonstructElmConc_rpvr => plt_biom%RootNonstructElmConc_rpvr, &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    icwood                    => pltpar%icwood,                      &
    iroot                     => pltpar%iroot,                       &
    inonstruct                => pltpar%inonstruct,                  &
    k_woody_litr              => pltpar%k_woody_litr,                &
    k_fine_litr               => pltpar%k_fine_litr,                 &
    CNRTS_pft                 => plt_allom%CNRTS_pft,                &
    CPRTS_pft                 => plt_allom%CPRTS_pft,                &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,    &
    RootBiomGrosYld_pft       => plt_allom%RootBiomGrosYld_pft,      &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,          &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,            &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,       &
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr,          &
    Root1stXNumL_pvr          => plt_morph%Root1stXNumL_pvr,         &
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft,          &
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft,        &
    NIXBotRootLayer_rpft      => plt_morph%NIXBotRootLayer_rpft,     &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,      &
    SoilBulkDensity_vr         => plt_soilchem%SoilBulkDensity_vr,     &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,    &
    C4PhotosynDowreg_brch     => plt_photo%C4PhotosynDowreg_brch     &
  )  
!   call SumRootBiome(NZ,mass_inital)

  IF(SoilBulkDensity_vr(L).GT.ZERO)THEN
    Root1stDepz2Surf=Root1stDepz_pft(N,NR,NZ)-CumSoilThickness_vr(0)
  ELSE
    Root1stDepz2Surf=Root1stDepz_pft(N,NR,NZ)
  ENDIF

  litrflxt=0._r8;RCO2flxt=0._r8
  !root pass the top surface of layer L, and the given root referenced by (N,NR) has not been updated
  !identify root tip
  IF(Root1stDepz2Surf.GT.CumSoilThickness_vr(L-1) .AND. iRootXsUpdateFlag(N,NR).EQ.ifalse)THEN
    Root1stXNumL_pvr(N,L,NZ)=Root1stXNumL_pvr(N,L,NZ)+RootPrimeAxsNum
    !primary roots
    IF(Root1stDepz2Surf.LE.CumSoilThickness_vr(L) .OR. L.EQ.MaxNumRootLays)THEN
      iRootXsUpdateFlag(N,NR)=itrue
!
!     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
!     ATTRIBUTED TO CURRENT AXIS
!
!     Root1stSink_pvr=primary root sink strength
!     RootSinkC_vr=total root sink strength
!     FRTN=fraction of primary root sink strength in axis
!
      IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
        FRTN=Root1stSink_pvr(N,L,NR)/RootSinkC_vr(N,L)
      ELSE
        FRTN=1.0_r8
      ENDIF
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
      CALL RootGroWaterDependence(N,L,NZ,WFNR,WFNRG)

!
!     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
      IF(RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
        CNPG=AMIN1(RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmn,N,L,NZ) &
          +RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CNKI),RootNonstructElmConc_rpvr(ielmp,N,L,NZ) &
          /(RootNonstructElmConc_rpvr(ielmp,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CPKI))
      ELSE
        CNPG=1.0_r8
      ENDIF
!
!     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     Rmaint1st_CO2=root maintenance respiration
!     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT1N=primary root N mass
!     TFN6_vr=temperature function for root maintenance respiration
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     fRootGrowPSISense=growth function of root water potential
!
! Rmaint1st_CO2 = AZMAX1(RmSpecPlant*RootMyco1stElm_raxs(ielmn,N,NR,NZ))*TFN6_vr(L)
      Rmaint1st_CO2 = AZMAX1(RmSpecPlant*RootMyco1stStrutElms_rpvr(ielmn,N,L,NR,NZ))*TFN6_vr(L)
      IF(is_root_shallow(iPlantRootProfile_pft(NZ)).OR.&
        iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
        Rmaint1st_CO2=Rmaint1st_CO2*fRootGrowPSISense
      ENDIF
!
!     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RNonstCO2_OUltd=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     
!     FRTN=fraction of primary root sink strength in axis
!     CPOOL=non-structural C mass
!     fTgrowRootP_vr=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!     fRootGrowPSISense=growth function of root water potential
!
      RNonstCO2_OUltd=AZMAX1(VMXC*FRTN*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)) &
        *fTgrowRootP_vr(L,NZ)*CNPG*C4PhotosynDowreg_brch(MainBranchNum_pft(NZ),NZ)*fRootGrowPSISense

      IF(Root1stDepz2Surf.GE.CumSoilThickness_vr(MaxNumRootLays))THEN
        RNonstCO2_OUltd=AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)
      ENDIF
!
!     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!     RCO2XMaint1st_OUltd,RCO2XMaint1st_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RNonstCO2_Oltd=RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
      RCO2XMaint1st_OUltd=RNonstCO2_OUltd-Rmaint1st_CO2
      RCO2XMaint1st_Oltd=RNonstCO2_Oltd-Rmaint1st_CO2
      RGrowCO2_OUltd=AZMAX1(RCO2XMaint1st_OUltd)*WFNRG
      RGrowCO2_Oltd=AZMAX1(RCO2XMaint1st_Oltd)*WFNRG
!
!     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS_pft,CPRTS_pft=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))
      PPOOLB=AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS_pft(NZ),PPOOLB*DMRTR/CPRTS_pft(NZ))
      IF(RGrowCO2_OUltd.GT.0.0_r8)THEN
        RGrowCO2_OUltd=AMIN1(RGrowCO2_OUltd,FNP)
      ELSE
        RGrowCO2_OUltd=0._r8
      ENDIF

      IF(RGrowCO2_Oltd.GT.0.0_r8)THEN
        RGrowCO2_Oltd=AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ))
      ELSE
        RGrowCO2_Oltd=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     RootMycoNonst4GrowC_OUltd,RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn unltd,ltd by O2
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     RootMycoNonst4Grow_OUltd(ielmc),RootMycoNonst4Grow_Oltd(ielmc)=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     RootMycoNonst4Grow_OUltd(ielmn),RootMycoNonst4Grow_Oltd(ielmn),RootMycoNonst4Grow_Oltd(ielmp)=nonstructural N,P unltd,ltd by O2 used in growth
!     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
!
      RootMycoNonst4GrowC_OUltd=RGrowCO2_OUltd/DMRTD
      RootMycoNonst4GrowC_Oltd=RGrowCO2_Oltd/DMRTD
      RootMycoNonst4Grow_OUltd(ielmc)=RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
      RootMycoNonst4Grow_OUltd(ielmn)=AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
      RootMycoNonst4Grow_OUltd(ielmp)=AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)
      RootMycoNonst4Grow_Oltd(ielmc)=RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)      
      RootMycoNonst4Grow_Oltd(ielmn)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmn,N,L,NZ) &
        ,RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
      RootMycoNonst4Grow_Oltd(ielmp)=AZMAX1(AMIN1(FRTN*RootMycoNonstElms_rpvr(ielmp,N,L,NZ) &
        ,RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))
      RCO2Nonst4Nassim_OUltd=AZMAX1(1.70_r8*RootMycoNonst4Grow_OUltd(ielmn))
      RCO2Nonst4Nassim_Oltd=AZMAX1(1.70_r8*RootMycoNonst4Grow_Oltd(ielmn))

!
!     TOTAL PRIMARY ROOT RESPIRATION
!
!     RCO2T1st_OUltd,RCO2T1st_Oltd=total C respiration unltd,ltd by O2
!     Rmaint1st_CO2=root maintenance respiration
!     RNonstCO2_OUltd,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
!     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
!     RootCO2Autor_pvr=total root respiration
!     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
!
      dRootMycoElms=0._r8
      RCO2T1st_OUltd=AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd+RCO2Nonst4Nassim_OUltd
      RCO2T1st_Oltd =AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+RCO2Nonst4Nassim_Oltd
      
      call RemobilizePrimeRoots(N,L,NZ,NR,RCO2XMaint1st_Oltd,RCO2XMaint1st_OUltd,dRootMycoElms,&
        RCO2T1st_OUltd,RCO2T1st_Oltd,Frac2Senes1,litrflxt)

      dRootMycoElms(ielmc)=dRootMycoElms(ielmc)-RCO2T1st_Oltd

      DO NE=1,NumPlantChemElms
        dRootMycoElms(NE)=dRootMycoElms(NE)-RootMycoNonst4Grow_Oltd(NE)
        RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+dRootMycoElms(NE)
      ENDDO         

!
!     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
!     THROUGH WHICH PRIMARY ROOTS GROW
!
!     RTDP1=primary root depth from soil surface
!     CumSoilThickness_vr=depth from soil surface to layer bottom
!     Root1stLen_rpvr=primary root length
!     SeedDepth_pft=seeding depth
!     FRCO2=fraction of primary root respiration attributed to layer
!     RootCO2Autor_pvr=total root respiration
!     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
!     RCO2T1st_OUltd,RCO2T1st_Oltd=total C respiration unltd,ltd by O2
!
      IF(Root1stDepz_pft(N,NR,NZ).GT.CumSoilThickness_vr(NGTopRootLayer_pft(NZ)))THEN
        TFRCO2=0._r8
        D5100: DO LL=NGTopRootLayer_pft(NZ),NIXBotRootLayer_rpft(NR,NZ)
          IF(LL.LT.NIXBotRootLayer_rpft(NR,NZ))THEN
            FRCO2=AMIN1(1.0_r8,Root1stLen_rpvr(N,LL,NR,NZ)/(Root1stDepz_pft(N,NR,NZ)-SeedDepth_pft(NZ)))
          ELSE
            FRCO2=1.0_r8-TFRCO2
          ENDIF
          TFRCO2=TFRCO2+FRCO2
          RootRespPotent_pvr(N,LL,NZ)=RootRespPotent_pvr(N,LL,NZ)+RCO2T1st_OUltd*FRCO2
          RootCO2EmisPot_pvr(N,LL,NZ)=RootCO2EmisPot_pvr(N,LL,NZ)+RCO2T1st_Oltd*FRCO2
          RootCO2Autor_pvr(N,LL,NZ)=RootCO2Autor_pvr(N,LL,NZ)-RCO2T1st_Oltd*FRCO2
          RCO2flxt=RCO2flxt-RCO2T1st_Oltd*FRCO2
        ENDDO D5100
      ELSE
        RootRespPotent_pvr(N,L,NZ)=RootRespPotent_pvr(N,L,NZ)+RCO2T1st_OUltd
        RootCO2EmisPot_pvr(N,L,NZ)=RootCO2EmisPot_pvr(N,L,NZ)+RCO2T1st_Oltd
        RootCO2Autor_pvr(N,L,NZ)=RootCO2Autor_pvr(N,L,NZ)-RCO2T1st_Oltd
        RCO2flxt=RCO2flxt-RCO2T1st_Oltd
      ENDIF
!
!     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
!     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
!     HAVE DISAPPEARED
!
!     RootMycoNonst4Grow_Oltd(ielmc)=primary root C growth ltd by O2
!     Root1stNetGrowthElms(:)=net primary root C,N,P growth
!     Frac2Senes1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RootMycoNonst4Grow_Oltd(ielmn),RootMycoNonst4Grow_Oltd(ielmp)=nonstructural N,P ltd by O2 used in growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     Root2ndLen_rpvr=secondary root length

      Root1stNetGrowthElms(:)=RootMycoNonst4Grow_Oltd(:)

      if(Frac2Senes1>0._r8)then
        !remobilization nonzero
        DO NE=1,NumPlantChemElms
          Root1stNetGrowthElms(NE)=Root1stNetGrowthElms(NE)-Frac2Senes1*RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)  
        ENDDO
      endif

      !negative growth
      IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN      
        call Withdraw2ndRoots(N,NZ,L,NR,Root1stNetGrowthElms,litrflxt)
      ENDIF

      call ExtendPrimeRoots(L,L1,N,NR,NZ,WFNR,FRTN,RootMycoNonst4Grow_Oltd(ielmc),Root1stNetGrowthElms,&
        TotRoot1stLen,Root1stC)
    ENDIF  

!    print*,'RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)',RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ),&
!      L,NIXBotRootLayer_rpft(NR,NZ)
    
    IF(L.EQ.NIXBotRootLayer_rpft(NR,NZ))THEN
 
      call WithdrawPrimeRoots(L,NR,NZ,N,Root1stDepz2Surf,RootSinkC_vr,Root1stSink_pvr &
        ,Root2ndSink_pvr,RootPrimeAxsNum)
    ENDIF

!     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
!   what if nonstructural C is insufficient to meet the negative root mass
    IF(RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
!      print*,'rm as nonst',RootMycoNonstElms_rpvr(ielmc,N,L,NZ),RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)

      RootMycoNonstElms_rpvr(ielmc,N,L,NZ)       = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
      RootMyco1stElm_raxs(ielmc,N,NR,NZ)         = RootMyco1stElm_raxs(ielmc,N,NR,NZ)+RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
      RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ) = 0._r8
    ENDIF

    IF(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
      RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ) &
        +RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
      RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)=0._r8
    ENDIF

    !correct negative value
    if(abs(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))<1.e-8_r8)RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=0._r8

!
!     TOTAL PRIMARY ROOT LENGTH AND MASS
!
!     Root1stLen_rpvr=primary root length in soil layer
!     WTRT1=primary root C mass in soil layer
!     NIXBotRootLayer_rpft=deepest root layer
!
!    TotRoot1stLen=TotRoot1stLen+Root1stLen_rpvr(N,L,NR,NZ)
!    Root1stC=Root1stC+RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
    NIXBotRootLayer_rpft(NR,NZ)=MIN(NIXBotRootLayer_rpft(NR,NZ),MaxNumRootLays)
    IF(L.EQ.NIXBotRootLayer_rpft(NR,NZ))NRX(N,NR)=itrue
  !   call SumRootBiome(NZ,mass_finale)  

  ENDIF
  end associate
  end subroutine Grow1stRootAxes
!------------------------------------------------------------------------------------------

  subroutine RootGroWaterDependence(N,L,NZ,WFNR,WFNRG)
  implicit none
  integer, intent(in) :: N,L,NZ
  real(r8), intent(out) :: WFNR
  real(r8), intent(out) :: WFNRG

  real(r8) :: SoilResit4PrimRootPentration  

  associate(                                                             &
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr,            &  
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,                  &    
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,        &    
    SoilResit4RootPentrate_vr => plt_soilchem%SoilResit4RootPentrate_vr  &
  )
!     SoilResit4RootPentrate_vr,SoilResit4PrimRootPentration=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
!     WFNR=water function for root extension
!     WFNRG=respiration function of root water potential
!
  SoilResit4PrimRootPentration = SoilResit4RootPentrate_vr(L)*Root1stRadius_pvr(N,L,NZ)/1.0E-03_r8
  WFNR                         = AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-PSIMin4OrganExtens-SoilResit4PrimRootPentration))
  WFNRG                        = fRespWatSens(WFNR,iPlantRootProfile_pft(NZ))

  end associate    
  END subroutine RootGroWaterDependence      
!------------------------------------------------------------------------------------------

  subroutine Withdraw2ndRoots(N,NZ,L,NR,RootNetGrowthElms,litrflxt)
  !
  !Description
  !Due to negative growth of primary roots, kill secondary roots to make up 
  !this deficit
  implicit none
  integer, intent(in) :: N,NZ,L,NR
  real(r8), intent(inout) :: RootNetGrowthElms(NumPlantChemElms)
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)

  integer :: LX,LL,NE,M
  real(r8) :: GRTWTM
  real(r8) :: FSNCM,FSNCP
  real(r8) :: dRootMyco2ndst2Litr(NumPlantChemElms)
  real(r8) :: dRootMycoNonst2Litr(NumPlantChemElms)
  real(r8) :: dsenecE,dwoodyE,dfineE1,dfineE2
  associate(                                                         &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr,   &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
    icwood                    => pltpar%icwood,                      &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,      &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,    &
    iroot                     => pltpar%iroot,                       &
    k_woody_litr              => pltpar%k_woody_litr,                &
    k_fine_litr               => pltpar%k_fine_litr,                 &
    inonstruct                => pltpar%inonstruct,                  &
    Root2ndLen_rpvr            => plt_morph%Root2ndLen_rpvr            &
  )  

  LX=MAX(1,L-1)        
  D5105: DO LL=L,LX,-1
    GRTWTM=RootNetGrowthElms(ielmc)
    !grow secondary roots
    D5106: DO NE=1,NumPlantChemElms
      IF(RootNetGrowthElms(NE).LT.0.0_r8)THEN
        IF(RootNetGrowthElms(NE).GT.-RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ))THEN
          if(NE==ielmc)then
            Root2ndLen_rpvr(N,LL,NR,NZ)=Root2ndLen_rpvr(N,LL,NR,NZ)&
              +RootNetGrowthElms(NE)*Root2ndLen_rpvr(N,LL,NR,NZ)/RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)
          endif   
          RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)&
            +RootNetGrowthElms(NE)
          RootNetGrowthElms(NE)=0._r8
        ELSE
          if(NE==ielmc)Root2ndLen_rpvr(N,LL,NR,NZ)=0._r8
          RootNetGrowthElms(NE)=RootNetGrowthElms(NE)+RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)
          RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)=0._r8
        ENDIF
      ENDIF
    ENDDO D5106
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
!     Root2ndLen_rpvr=mycorrhizal length
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
!
    IF(GRTWTM.LT.0.0_r8)THEN
      IF(RootMyco2ndStrutElms_rpvr(ielmc,ipltroot,LL,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FSNCM=AMIN1(1.0_r8,ABS(GRTWTM)/RootMyco2ndStrutElms_rpvr(ielmc,ipltroot,LL,NR,NZ))
      ELSE
        FSNCM=1.0_r8
      ENDIF

      IF(RootMycoActiveBiomC_pvr(ipltroot,LL,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FSNCP=AMIN1(1.0_r8,ABS(GRTWTM)/RootMycoActiveBiomC_pvr(ipltroot,LL,NZ))
      ELSE
        FSNCP=1.0_r8
      ENDIF

      dRootMyco2ndst2Litr=0._r8
      dRootMycoNonst2Litr=0._r8
      D6450: DO M=1,jsken
        D64511: DO NE=1,NumPlantChemElms
          dsenecE=FSNCM*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ))
          dwoodyE=ElmAllocmat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_woody_litr)
          dfineE1=ElmAllocmat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_fine_litr)
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,LL,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,LL,NZ)+dwoodyE
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ)+dfineE1
          dRootMyco2ndst2Litr(NE)=dRootMyco2ndst2Litr(NE)+dwoodyE+dfineE1

          dfineE2=ElmAllocmat4Litr(NE,inonstruct,M,NZ)*FSNCP*(RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ))
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,LL,NZ)+dfineE2
          dRootMycoNonst2Litr(NE)=dRootMycoNonst2Litr(NE)+dfineE2

          litrflxt(NE)=litrflxt(NE)+dwoodyE+dfineE1+dfineE2
                            
        ENDDO D64511   
      ENDDO D6450

      DO NE=1,NumPlantChemElms  
        RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ) &
          -dRootMyco2ndst2Litr(NE)

        RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ)=RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ) &
          -dRootMycoNonst2Litr(NE)
        
      ENDDO
      Root2ndLen_rpvr(imycorrhz,LL,NR,NZ)=AZMAX1(Root2ndLen_rpvr(imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)
    ENDIF
  ENDDO D5105
  end associate
  end subroutine Withdraw2ndRoots        

!------------------------------------------------------------------------------------------
  subroutine RemobilizePrimeRoots(N,L,NZ,NR,RCO2XMaint1st_Oltd,RCO2XMaint1st_OUltd,dRootMycoElms,&
    RCO2T1st_OUltd,RCO2T1st_Oltd,Frac2Senes1,litrflxt)

  implicit none
  integer , intent(in) :: N,L,NZ,NR
  real(r8), intent(in) :: RCO2XMaint1st_Oltd       !diagnostic for maintenance deficit
  real(r8), intent(in) :: RCO2XMaint1st_OUltd  
  real(r8), intent(inout) :: RCO2T1st_Oltd
  real(r8), intent(inout) :: RCO2T1st_OUltd
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: dRootMycoElms(NumPlantChemElms)
  real(r8), intent(out) :: Frac2Senes1    
  integer :: M,NE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: Root1stStrutRemob(NumPlantChemElms)

  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: RCO2Nonst4Xmaint1st_Oltd,RCO2Nonst4Xmaint1st_OUltd
  real(r8) :: dsenecE
! begin_execution
  associate(                                                         &
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,       &
    RootNonstructElmConc_rpvr  => plt_biom%RootNonstructElmConc_rpvr,  &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    ZERO                      => plt_site%ZERO,                      &
    icwood                    => pltpar%icwood,                      &
    k_woody_litr              => pltpar%k_woody_litr,                &
    k_fine_litr               => pltpar%k_fine_litr,                 &
    iroot                     => pltpar%iroot,                       &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,    &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,      &
    RAutoRootO2Limter_rpvr     => plt_rbgc%RAutoRootO2Limter_rpvr,     &
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch,      &
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft         &
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
    .AND.RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmn,N,L,NZ)&
      ,RootNonstructElmConc_rpvr(ielmn,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CNKI) &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmp,N,L,NZ),RootNonstructElmConc_rpvr(ielmp,N,L,NZ)&
        +RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ),&
      RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ),&
      RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/CPKI)))
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
!     RCO2XMaint1st_OUltd,RCO2XMaint1st_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2Nonst4Xmaint1st_OUltd,RCO2Nonst4Xmaint1st_Oltd=excess maintenance respiration unltd,ltd by O2
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!     Root1stStrutRemob(ielmc),Root1stStrutRemob(ielmn),Root1stStrutRemob(ielmp)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     Frac2Senes1=fraction of primary root C to be remobilized
!
  IF(-RCO2XMaint1st_OUltd.GT.0.0_r8)THEN
    IF(-RCO2XMaint1st_OUltd.LT.RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
      RCO2Nonst4Xmaint1st_OUltd=-RCO2XMaint1st_OUltd
    ELSE
      RCO2Nonst4Xmaint1st_OUltd=AZMAX1(RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)
    ENDIF
  ELSE
    RCO2Nonst4Xmaint1st_OUltd=0._r8
  ENDIF

  IF(-RCO2XMaint1st_Oltd.GT.0.0_r8)THEN
    !maintenance deficit is less than remobilizable C.
    IF(-RCO2XMaint1st_Oltd.LT.RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
      RCO2Nonst4Xmaint1st_Oltd=-RCO2XMaint1st_Oltd
    ELSE
      RCO2Nonst4Xmaint1st_Oltd=AZMAX1(RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)*RAutoRootO2Limter_rpvr(N,L,NZ)
    ENDIF
  ELSE
    RCO2Nonst4Xmaint1st_Oltd=0._r8
  ENDIF

  IF(RCO2Nonst4Xmaint1st_Oltd.GT.0.0_r8 .AND. RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    Root1stStrutRemob(ielmc) = RCCC*RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
    Root1stStrutRemob(ielmn) = RootMyco1stStrutElms_rpvr(ielmn,N,L,NR,NZ)*(RCCN+(1.0_r8-RCCN) &
      *Root1stStrutRemob(ielmc)/RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ))      
    Root1stStrutRemob(ielmp)=RootMyco1stStrutElms_rpvr(ielmp,N,L,NR,NZ)*(RCCP+(1.0_r8-RCCP) &
      *Root1stStrutRemob(ielmc)/RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ))

    IF(Root1stStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
      Frac2Senes1=AZMAX1(AMIN1(1.0_r8,RCO2Nonst4Xmaint1st_Oltd/Root1stStrutRemob(ielmc)))
    ELSE
      Frac2Senes1=1.0_r8
    ENDIF
  ELSE
    Root1stStrutRemob(1:NumPlantChemElms)=0._r8
    Frac2Senes1=0._r8
  ENDIF
!
!     PRIMARY ROOT LitrFall CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     Frac2Senes1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     Root1stStrutRemob(ielmc),Root1stStrutRemob(ielmn),Root1stStrutRemob(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
  if(Frac2Senes1>0._r8)then
    D6355: DO M=1,jsken
      DO NE=1,NumPlantChemElms    
        dsenecE                                      = Frac2Senes1*AZMAX1(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)-Root1stStrutRemob(NE))
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ) = LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ) &
          +ElmAllocmat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_woody_litr)
        litrflxt(NE)=litrflxt(NE)+ElmAllocmat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_woody_litr)
        
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ) &
          +ElmAllocmat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_fine_litr)
        litrflxt(NE)=litrflxt(NE)+ElmAllocmat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAlloc2Litr(NE,k_fine_litr)  
      ENDDO
    ENDDO D6355
  endif
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     Rmaint1st_CO2=root maintenance respiration
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn ltd by O2
!     RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
!     RCO2Nonst4Xmaint1st_Oltd=excess maintenance respiration ltd by O2
!     Frac2Senes1=fraction of primary root C to be remobilized
!     Root2ndStrutRemob(:)=remobilization of C,N,P from senescing root
!     RootMycoNonst4Grow_Oltd(ielmn),RootMycoNonst4Grow_Oltd(ielmp)=nonstructural N,P ltd by O2 used in growth
!  
  RCO2T1st_Oltd=RCO2T1st_Oltd+RCO2Nonst4Xmaint1st_Oltd
  RCO2T1st_OUltd=RCO2T1st_OUltd+RCO2Nonst4Xmaint1st_OUltd
  DO NE=1,NumPlantChemElms
    dRootMycoElms(NE)=dRootMycoElms(NE)+Frac2Senes1*Root1stStrutRemob(NE)
  ENDDO 
  end associate
  end subroutine RemobilizePrimeRoots

!------------------------------------------------------------------------------------------
  subroutine ExtendPrimeRoots(L,L1,N,NR,NZ,WFNR,FRTN,RootMycoNonstC4Grow_Oltd,RootNetGrowthElms,&
    TotRoot1stLen,Root1stC)
  !
  !
  implicit none
  integer, intent(in) :: L,L1
  integer, intent(in) :: N    !root/myco indicator
  integer, intent(in) :: NR   !root axis 
  integer, intent(in) :: NZ   !pft number
  real(r8), intent(in):: WFNR,FRTN
  real(r8), intent(in) :: RootMycoNonstC4Grow_Oltd !oxygen limited root C yield for growth
  real(r8), intent(in) :: RootNetGrowthElms(NumPlantChemElms)
  real(r8), intent(inout) :: TotRoot1stLen,Root1stC
  real(r8) :: Root1stExtension
  real(r8) :: FGROL,FGROZ
  integer :: NE
  REAL(R8) :: XFRE(NumPlantChemElms)
! begin_execution
  associate(                                                         &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,       &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,          &
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,        &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    rCNNonstRemob_pft         => plt_allom%rCNNonstRemob_pft,        &
    rCPNonstRemob_pft         => plt_allom%rCPNonstRemob_pft,        &
    FracRootElmAlloc2Litr     => plt_allom%FracRootElmAlloc2Litr,    &
    CumSoilThickness_vr       => plt_site%CumSoilThickness_vr,       &
    DLYR3                     => plt_site%DLYR3,                     &
    MaxNumRootLays            => plt_site%MaxNumRootLays,            &
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,       &
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,              &
    PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr,              &
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                 &
    k_woody_litr              => pltpar%k_woody_litr,                &
    k_fine_litr               => pltpar%k_fine_litr,                 &
    Root1stSpecLen_pft        => plt_morph%Root1stSpecLen_pft,       &
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr,        &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,          &
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft,          &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,       &
    NIXBotRootLayer_rpft      => plt_morph%NIXBotRootLayer_rpft,     &
    SeedDepth_pft             => plt_morph%SeedDepth_pft             &
  )
!     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     Root1stExtension=primary root length extension
!     RootMycoNonst4Grow_Oltd=primary root C growth ltd by O2
!     Root1stSpecLen_pft=specific primary root length from startq.f
!     PP=PFT population
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     RootNetGrowthElms(:)=net primary root C,N,P growth
!     RTDP1=primary root depth from soil surface
!     SeedDepth_pft=seeding depth
!     Frac2Senes1=fraction of primary root C to be remobilized
!     Root1stLen_rpvr=primary root length
!     RootNetGrowthElms(:)=net root C,N,P growth
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     DLYR=soil layer thickness
!
  Root1stExtension=RootMycoNonstC4Grow_Oltd*Root1stSpecLen_pft(N,NZ)/PlantPopulation_pft(NZ)&
    *WFNR*FracRootElmAlloc2Litr(ielmc,k_fine_litr)
  IF(RootNetGrowthElms(ielmc).LT.0.0_r8 .AND. RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    !primary roots withdraw
    Root1stExtension=Root1stExtension+RootNetGrowthElms(ielmc) &
      *(Root1stDepz_pft(N,NR,NZ)-SeedDepth_pft(NZ))/RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
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
! Question, 03/29/2024, Jinyun Tang: needs double check the following calculation 
! if FGROL < 1.0, then the extension is all in current layer, meaning FGROZ=0.0
!  
  IF(Root1stExtension.GT.ZERO4Groth_pft(NZ) .AND. L.LT.MaxNumRootLays)THEN
    FGROL=AZMAX1(AMIN1(1.0_r8,(CumSoilThickness_vr(L)-Root1stDepz_pft(N,NR,NZ))/Root1stExtension))
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
!     RootNetGrowthElms(:)=net root C,N,P growth
!     Root1stExtension=primary root length extension
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     FGROL,FGROZ=fraction of Root1stExtension in current,next lower soil layer
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
!     Root1stLen_rpvr=primary root length
!
  Root1stDepz_pft(N,NR,NZ)=Root1stDepz_pft(N,NR,NZ)+Root1stExtension

  DO NE=1,NumPlantChemElms
    RootMyco1stElm_raxs(NE,N,NR,NZ)         = RootMyco1stElm_raxs(NE,N,NR,NZ)+RootNetGrowthElms(NE)
    RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)+RootNetGrowthElms(NE)*FGROL 
  ENDDO

  RootProteinC_pvr(N,L,NZ) = RootProteinC_pvr(N,L,NZ)+ &
    AMIN1(rCNNonstRemob_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmn,N,L,NR,NZ) &
         ,rCPNonstRemob_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmp,N,L,NR,NZ))
  Root1stLen_rpvr(N,L,NR,NZ)=Root1stLen_rpvr(N,L,NR,NZ)+Root1stExtension*FGROL
!
!     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
!     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
!     OF CURRENT SOIL LAYER
!
!     FGROZ=fraction of Root1stExtension in next lower soil layer
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     RootNetGrowthElms(:)=net root C,N,P growth
!     RootProteinC_pvr=total root protein C mass
!     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
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
      RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ)+RootNetGrowthElms(NE)*FGROZ
    ENDDO

    RootProteinC_pvr(N,L1,NZ) = RootProteinC_pvr(N,L1,NZ) &
      +AMIN1(rCNNonstRemob_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmn,N,L1,NR,NZ) &
            ,rCPNonstRemob_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmp,N,L1,NR,NZ))
    PopuRootMycoC_pvr(N,L1,NZ)  = PopuRootMycoC_pvr(N,L1,NZ)+RootMyco1stStrutElms_rpvr(ielmc,N,L1,NR,NZ)
    Root1stLen_rpvr(N,L1,NR,NZ) = Root1stLen_rpvr(N,L1,NR,NZ)+Root1stExtension*FGROZ
    Root1stRadius_pvr(N,L1,NZ)  = Root1stRadius_pvr(N,L,NZ)
    TotRoot1stLen               = TotRoot1stLen+Root1stLen_rpvr(N,L1,NR,NZ)
    Root1stC                    = Root1stC+RootMyco1stStrutElms_rpvr(ielmc,N,L1,NR,NZ)

    DO NE=1,NumPlantChemElms
      XFRE(NE)                           = FRTN* RootMycoNonstElms_rpvr(NE,N,L,NZ)
      RootMycoNonstElms_rpvr(NE,N,L,NZ)  = RootMycoNonstElms_rpvr(NE,N,L,NZ)-XFRE(NE)
      RootMycoNonstElms_rpvr(NE,N,L1,NZ) = RootMycoNonstElms_rpvr(NE,N,L1,NZ)+XFRE(NE)
    ENDDO
    PSIRoot_pvr(N,L1,NZ)        = PSIRoot_pvr(N,L,NZ)
    PSIRootOSMO_vr(N,L1,NZ)     = PSIRootOSMO_vr(N,L,NZ)
    PSIRootTurg_vr(N,L1,NZ)     = PSIRootTurg_vr(N,L,NZ)
    NIXBotRootLayer_rpft(NR,NZ) = MAX(NGTopRootLayer_pft(NZ),L+1)
  ENDIF
  end associate
  end subroutine ExtendPrimeRoots
!------------------------------------------------------------------------------------------

  subroutine WithdrawPrimeRoots(L,NR,NZ,N,Root1stDepz2Surf,RootSinkC_vr &
    ,Root1stSink_pvr,Root2ndSink_pvr,RootPrimeAxsNum)
  implicit none
  integer, intent(in) :: L,NR,NZ,N
  real(r8), intent(in):: Root1stDepz2Surf
  real(r8), INTENT(IN) :: RootSinkC_vr(jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: RootPrimeAxsNum
  integer :: LL,NN,NE,NTG
  real(r8) :: XFRD,XFRW,FRTN
  real(r8) :: XFRE
  logical :: RootDepzChk

! begin_execution
  associate(                                                         &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr, &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,    &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,          &
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,        &
    RootNodulStrutElms_rpvr    => plt_biom%RootNodulStrutElms_rpvr,    &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    RootNodulNonstElms_rpvr    => plt_biom%RootNodulNonstElms_rpvr,    &
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft,    &
    CumSoilThickness_vr       => plt_site%CumSoilThickness_vr,       &
    ZEROS2                    => plt_site%ZEROS2,                    &
    DLYR3                     => plt_site%DLYR3,                     &
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr,           &
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr,           &
    VLSoilPoreMicP_vr         => plt_soilchem%VLSoilPoreMicP_vr,     &
    Root1stXNumL_pvr          => plt_morph%Root1stXNumL_pvr,         &
    Root2ndXNum_pvr           => plt_morph%Root2ndXNum_pvr,          &
    Root1stLen_rpvr           => plt_morph%Root1stLen_rpvr,          &
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr,         &
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft,          &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,       &
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,       &
    MY                        => plt_morph%MY,                       &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,            &
    NIXBotRootLayer_rpft      => plt_morph%NIXBotRootLayer_rpft      &
  )
!     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
!     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
!     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
!     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
!
!     NIXBotRootLayer_rpft=deepest root layer
!     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
!     Root1stDepz2Surf=primary root depth from soil surface
!     CumSoilThickness_vr=depth from soil surface to layer bottom
!     SeedDepth_pft=seeding depth
!     FRTN=fraction of primary root sink strength in axis
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     Root1stLen_rpvr=primary root length
!     RootProteinC_pvr=root protein C mass
!     WTRTD=root C mass
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root

  D5115: DO LL=L,NGTopRootLayer_pft(NZ)+1,-1
    RootDepzChk=Root1stDepz2Surf.LT.CumSoilThickness_vr(LL-1) .OR. Root1stDepz2Surf.LT.SeedDepth_pft(NZ)
    IF(VLSoilPoreMicP_vr(LL-1).GT.ZEROS2 .AND. RootDepzChk)THEN
      IF(RootSinkC_vr(N,LL).GT.ZERO4Groth_pft(NZ))THEN
        FRTN=(Root1stSink_pvr(N,LL,NR)+Root2ndSink_pvr(N,LL,NR))/RootSinkC_vr(N,LL)
      ELSE
        FRTN=1.0_r8
      ENDIF
!      print*,'begd5110'
      D5110: DO NN=1,MY(NZ)
        Root1stLen_rpvr(NN,LL-1,NR,NZ)=Root1stLen_rpvr(NN,LL-1,NR,NZ)+Root1stLen_rpvr(NN,LL,NR,NZ)
        Root1stLen_rpvr(NN,LL,NR,NZ)  =0._r8
        DO NE=1,NumPlantChemElms
          RootMyco1stStrutElms_rpvr(NE,NN,LL-1,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,NN,LL-1,NR,NZ) &
            +RootMyco1stStrutElms_rpvr(NE,NN,LL,NR,NZ)
          RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ) &
            +RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ)

          RootMyco1stStrutElms_rpvr(NE,NN,LL,NR,NZ) = 0._r8
          RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ) = 0._r8
          XFRE                                      = FRTN*RootMycoNonstElms_rpvr(NE,NN,LL,NZ)
          RootMycoNonstElms_rpvr(NE,NN,LL,NZ)       = RootMycoNonstElms_rpvr(NE,NN,LL,NZ)-XFRE
          RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ)     = RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ)+XFRE
        ENDDO
        XFRW = FRTN*RootProteinC_pvr(NN,L,NZ)
        XFRD = FRTN* PopuRootMycoC_pvr(NN,LL,NZ)

        RootProteinC_pvr(NN,LL,NZ)    = RootProteinC_pvr(NN,LL,NZ)-XFRW
        RootProteinC_pvr(NN,LL-1,NZ)  = RootProteinC_pvr(NN,LL-1,NZ)+XFRW
        PopuRootMycoC_pvr(NN,LL,NZ)   = PopuRootMycoC_pvr(NN,LL,NZ)-XFRD
        PopuRootMycoC_pvr(NN,LL-1,NZ) = PopuRootMycoC_pvr(NN,LL-1,NZ)+XFRD
!
!     WITHDRAW GASES IN PRIMARY ROOTS
!
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
!     FRTN=fraction of primary root sink strength in axis
!
        DO NTG=idg_beg,idg_NH3
          RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-FRTN &
            *(trcg_rootml_pvr(idg_CO2,NN,LL,NZ)+trcs_rootml_pvr(idg_CO2,NN,LL,NZ))
          trcg_rootml_pvr(NTG,NN,LL,NZ) = (1.0_r8-FRTN)*trcg_rootml_pvr(NTG,NN,LL,NZ)
          trcs_rootml_pvr(NTG,NN,LL,NZ) = (1.0_r8-FRTN)*trcs_rootml_pvr(NTG,NN,LL,NZ)
        ENDDO
      ENDDO D5110
!
!     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
!
!     Root2ndXNum_rpvr,Root2ndXNum_pvr=number of secondary root axes
!     Root1stXNumL_pvr=number of primary root axes
!     Root1stLen_rpvr=primary root length
!     CumSoilThickness_vr=depth from soil surface to layer bottom
!     SeedDepth_pft=seeding depth
!
      Root2ndXNum_pvr(N,LL,NZ)=Root2ndXNum_pvr(N,LL,NZ)-Root2ndXNum_rpvr(N,LL,NR,NZ)
      Root2ndXNum_pvr(N,LL-1,NZ)=Root2ndXNum_pvr(N,LL-1,NZ)+Root2ndXNum_rpvr(N,LL,NR,NZ)
      Root2ndXNum_rpvr(N,LL,NR,NZ)=0._r8
      Root1stXNumL_pvr(N,LL,NZ)=Root1stXNumL_pvr(N,LL,NZ)-RootPrimeAxsNum
      IF(LL-1.GT.NGTopRootLayer_pft(NZ))THEN
        Root1stLen_rpvr(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness_vr(LL-1)-Root1stDepz_pft(N,NR,NZ))
      ELSE
        Root1stLen_rpvr(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness_vr(LL-1)-Root1stDepz_pft(N,NR,NZ)) &
          -(SeedDepth_pft(NZ)-CumSoilThickness_vr(LL-2))
      ENDIF
!
!     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
!
!     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     FRTN=fraction of primary root sink strength in axis
!     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
!
      IF(is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
        DO NE=1,NumPlantChemElms
          XFRE=FRTN*RootNodulStrutElms_rpvr(NE,LL,NZ)
          RootNodulStrutElms_rpvr(NE,LL,NZ)=RootNodulStrutElms_rpvr(NE,LL,NZ)-XFRE
          RootNodulStrutElms_rpvr(NE,LL-1,NZ)=RootNodulStrutElms_rpvr(NE,LL-1,NZ)+XFRE
          XFRE=FRTN*RootNodulNonstElms_rpvr(NE,LL,NZ)
          RootNodulNonstElms_rpvr(NE,LL,NZ)=RootNodulNonstElms_rpvr(NE,LL,NZ)-XFRE
          RootNodulNonstElms_rpvr(NE,LL-1,NZ)=RootNodulNonstElms_rpvr(NE,LL-1,NZ)+XFRE
        ENDDO
      ENDIF
      NIXBotRootLayer_rpft(NR,NZ)=MAX(NGTopRootLayer_pft(NZ),LL-1)
    ELSE
      EXIT
    ENDIF
  ENDDO D5115
  end associate
  end subroutine WithdrawPrimeRoots

!------------------------------------------------------------------------------------------
  subroutine SummarizeRootSink(I,J,NZ,RootPrimeAxsNum,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in):: RootPrimeAxsNum
  real(r8),INTENT(OUT) :: RootSinkC_vr(jroots,JZ1)
  real(r8),intent(out) :: Root1stSink_pvr(jroots,JZ1,NumOfCanopyLayers1)
  real(r8),intent(out) :: Root2ndSink_pvr(jroots,JZ1,NumOfCanopyLayers1)
  real(r8),INTENT(OUT) :: RootSinkC(jroots)
  integer :: N,L,K,NR,NE,ntu
  REAL(R8) :: Root1stLocDepz_vr(NumOfCanopyLayers1,JZ1)
  real(r8) :: RecoRootMycoC4Nup,CUPRO,CUPRC
  real(r8) :: RTDPP,RTDPS,RTSKP
  real(r8) :: RTSKS,rscal
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: RCO2flx

  associate(                                                   &
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr, &
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft,         &
    iPlantRootProfile_pft  => plt_pheno%iPlantRootProfile_pft, &
    VLSoilPoreMicP_vr      => plt_soilchem%VLSoilPoreMicP_vr,  &
    ZEROS2                 => plt_site%ZEROS2,                 &
    NU                     => plt_site%NU,                     &
    CumSoilThickness_vr    => plt_site%CumSoilThickness_vr,    &
    ZERO                   => plt_site%ZERO,                   &
    DLYR3                  => plt_site%DLYR3,                  &
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr,      &
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr,  &
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr,  &
    RootMycoExudEUptk_pvr    => plt_rbgc%RootMycoExudEUptk_pvr,    &
    RootCO2EmisPot_pvr     => plt_rbgc%RootCO2EmisPot_pvr,     &
    RootRespPotent_pvr     => plt_rbgc%RootRespPotent_pvr,     &
    RootCO2Autor_pvr       => plt_rbgc%RootCO2Autor_pvr,       &
    CanopyHeight4WatUptake_pft   => plt_morph%CanopyHeight4WatUptake_pft,  &
    MY                     => plt_morph%MY,                    &
    Root1stRadius_pvr      => plt_morph%Root1stRadius_pvr,     &
    Root1stDepz_pft        => plt_morph%Root1stDepz_pft,       &
    HypoctoHeight_pft      => plt_morph%HypoctoHeight_pft,     &
    Root2ndRadius_pvr      => plt_morph%Root2ndRadius_pvr,     &
    Root2ndXNum_rpvr       => plt_morph%Root2ndXNum_rpvr,      &
    Root2ndAveLen_pvr      => plt_morph%Root2ndAveLen_pvr,     &
    SeedDepth_pft          => plt_morph%SeedDepth_pft,         &
    MaxSoiL4Root_pft       => plt_morph%MaxSoiL4Root_pft,      &
    NumRootAxes_pft        => plt_morph%NumRootAxes_pft        &
  )

!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER

  RootSinkC_vr    = 0._R8
  Root1stSink_pvr = 0._r8
  Root2ndSink_pvr = 0._r8
  RootSinkC       = 0._r8
  RCO2flx         = 0._r8
!   call SumRootBiome(NZ,mass_inital)

  D4995: DO N=1,MY(NZ)
    D4990: DO L=NU,MaxSoiL4Root_pft(NZ)
!
!     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
!     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
!
!     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
!     RecoRootMycoC4Nup=C respiration for nutrient uptake
!     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
!     RootNH4Uptake_pvr,RootNH4BUptake_pvr,RUPN03,RootNO3BUptake_pvr=uptake from non-band,band of NH4,NO3
!     RootH2PO4Uptake_pvr,RootH2PO4BUptake_pvr,RootHPO4Uptake_pvr,RootNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4
!     RootOUlmNutUptake_pvr,RootOUlmNutUptake_pvr,RUON03,RootOUlmNutUptake_pvr=uptake from non-band,band of NH4,NO3 unlimited by O2
!     RootOUlmNutUptake_pvr,RootOUlmNutUptake_pvr,RootOUlmNutUptake_pvr,RootOUlmNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4 unlimited by O2
!     RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr,RUCN03,RootCUlmNutUptake_pvr=uptake from non-band,band of NH4,NO3 unlimited by nonstructural C
!     RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr,RootCUlmNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4 unlimited by nonstructural C
!     why is 0.86? it refers to C cost for N asimilation
!      print*,'xVLSoilPoreMicP_vr',L,MaxSoiL4Root_pft(NZ),VLSoilPoreMicP_vr(L),RecoRootMycoC4Nup
      IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
        RecoRootMycoC4Nup=0._r8;CUPRO=0._r8;CUPRC=0._r8
        do ntu=ids_nutb_beg+1,ids_nuts_end
          if(ntu==ids_NO2B .or. ntu==ids_NO2)cycle
          RecoRootMycoC4Nup = RecoRootMycoC4Nup+RootNutUptake_pvr(ntu,N,L,NZ)
          CUPRO             = CUPRO+RootOUlmNutUptake_pvr(ntu,N,L,NZ)
          CUPRC             = CUPRC+RootCUlmNutUptake_pvr(ntu,N,L,NZ)
        enddo  
        RecoRootMycoC4Nup = 0.86_r8*RecoRootMycoC4Nup
        CUPRO             = 0.86_r8*CUPRO
        CUPRC             = 0.86_r8*CUPRC
!
!       print*,'ACCUMULATE RESPIRATION IN FLUX ARRAYS'
!
!     RootCO2Autor_pvr=total root respiration
!     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
!     RecoRootMycoC4Nup=C respiration for nutrient uptake
!     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
!     CPOOLR=non-structural C mass in root
!
        RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+CUPRO
        RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+CUPRC

        if(abs(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))<1.e-10_r8)RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=0._r8
        if(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)<RecoRootMycoC4Nup .and. RecoRootMycoC4Nup>0._r8)then

          rscal=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)/RecoRootMycoC4Nup*0.99999_r8          
          
          do ntu=ids_nutb_beg+1,ids_nuts_end
            if(ntu==ids_NO2B .or. ntu==ids_NO2)cycle
            RootNutUptake_pvr(ntu,N,L,NZ) = RootNutUptake_pvr(ntu,N,L,NZ)*rscal
          enddo
          RecoRootMycoC4Nup=RecoRootMycoC4Nup*rscal
        endif

        RootCO2Autor_pvr(N,L,NZ)             = RootCO2Autor_pvr(N,L,NZ)-RecoRootMycoC4Nup
        RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RecoRootMycoC4Nup
        RCO2flx                              = RCO2flx-RecoRootMycoC4Nup
!
!     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
!     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
!     RootNH4Uptake_pvr,RootNH4BUptake_pvr,RUPN03,RootNO3BUptake_pvr=uptake from non-band,band of NH4,NO3
!     RootH2PO4Uptake_pvr,RootH2PO4BUptake_pvr,RootHPO4Uptake_pvr,RootNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4
!

        RootMycoNonstElms_rpvr(ielmn,N,L,NZ)=RootMycoNonstElms_rpvr(ielmn,N,L,NZ)+&
          (RootNutUptake_pvr(ids_NH4,N,L,NZ) +RootNutUptake_pvr(ids_NH4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_NO3,N,L,NZ)+RootNutUptake_pvr(ids_NO3B,N,L,NZ))
         RootMycoNonstElms_rpvr(ielmp,N,L,NZ)=RootMycoNonstElms_rpvr(ielmp,N,L,NZ) &
           +(RootNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootNutUptake_pvr(ids_H2PO4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ))

!
!       print*,'GROWTH OF EACH ROOT AXIS',N,L,MaxSoiL4Root_pft(NZ)
!
        D4985: DO NR=1,NumRootAxes_pft(NZ)
!
!     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
!
!     RTDP1=primary root depth from soil surface
!     RTDPP=primary root depth from canopy
!     CumSoilThickness_vr=depth from soil surface to layer bottom
!     CanopyHeight4WatUptake_pft=canopy height for water uptake
!     RTSK=relative primary root sink strength
!     Root1stSink_pvr=primary root sink strength
!     RootPrimeAxsNum=number of primary root axes
!     RRAD1,Root2ndRadius_pvr=primary, secondary root radius
!     RootSinkC,RootSinkC_vr=total root sink strength
!
          IF(N.EQ.ipltroot)THEN
            IF(Root1stDepz_pft(N,NR,NZ).GT.CumSoilThickness_vr(L-1))THEN
              IF(Root1stDepz_pft(N,NR,NZ).LE.CumSoilThickness_vr(L))THEN
                RTDPP                   = Root1stDepz_pft(N,NR,NZ)+CanopyHeight4WatUptake_pft(NZ)
                Root1stSink_pvr(N,L,NR) = RTSK(iPlantRootProfile_pft(NZ))&
                  *RootPrimeAxsNum*Root1stRadius_pvr(N,L,NZ)**2._r8/RTDPP
                RootSinkC(N)      = RootSinkC(N)+Root1stSink_pvr(N,L,NR)
                RootSinkC_vr(N,L) = RootSinkC_vr(N,L)+Root1stSink_pvr(N,L,NR)
              ENDIF
            ENDIF

!
!     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
!     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
!
!     Root1stLocDepz_vr=depth of primary root axis in layer
!     RTDP1=primary root depth from soil surface
!     CumSoilThickness_vr=depth from soil surface to layer bottom
!     RTDPX=distance behind growing point for secondary roots
!     DLYR=layer thickness
!     SeedDepth_pft=seeding depth
!     HypoctoHeight_pft=hypocotyledon height
!     CanopyHeight4WatUptake_pft=canopy height for water uptake
!     RTDPS=secondary root depth from canopy
!     RTSKP,RTSKS=primary,secondary root sink strength
!     RTN2=number of secondary root axes
!     Root2ndSink_pvr=total secondary root sink strength
!     Root2ndAveLen_pvr=average secondary root length
!     RootSinkC,RootSinkC_vr=total root sink strength
!
            Root1stLocDepz_vr(NR,L) = AZMAX1(Root1stDepz_pft(ipltroot,NR,NZ)-CumSoilThickness_vr(L-1)-RTDPX)
            Root1stLocDepz_vr(NR,L) = AZMAX1(AMIN1(DLYR3(L),Root1stLocDepz_vr(NR,L)) &
              -AZMAX1(SeedDepth_pft(NZ)-CumSoilThickness_vr(L-1)-HypoctoHeight_pft(NZ)))
            RTDPS=AMAX1(SeedDepth_pft(NZ),CumSoilThickness_vr(L-1))+0.5_r8*Root1stLocDepz_vr(NR,L)+CanopyHeight4WatUptake_pft(NZ)

            IF(RTDPS.GT.ZERO)THEN
              RTSKP = RootPrimeAxsNum*Root1stRadius_pvr(N,L,NZ)**2._r8/RTDPS
              RTSKS = safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)*Root2ndRadius_pvr(N,L,NZ)**2._r8,Root2ndAveLen_pvr(N,L,NZ))
              IF(RTSKP+RTSKS.GT.ZERO4Groth_pft(NZ))THEN
                Root2ndSink_pvr(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
              ELSE
                Root2ndSink_pvr(N,L,NR)=0._r8
              ENDIF
            ELSE
              Root2ndSink_pvr(N,L,NR)=0._r8
            ENDIF
          ELSE
            Root2ndSink_pvr(N,L,NR)=safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)&
              *Root2ndRadius_pvr(N,L,NZ)**2._r8,Root2ndAveLen_pvr(N,L,NZ))
          ENDIF

          RootSinkC(N)      = RootSinkC(N)+Root2ndSink_pvr(N,L,NR)
          RootSinkC_vr(N,L) = RootSinkC_vr(N,L)+Root2ndSink_pvr(N,L,NR)

        ENDDO D4985
      ENDIF
    ENDDO D4990
  ENDDO D4995

!   call SumRootBiome(NZ,mass_finale)

  end associate
  end subroutine SummarizeRootSink
!------------------------------------------------------------------------------------------
  subroutine RootCheck(I,J,NZ,info)
  !
  !root axes are sorted by root order and then also by vertical layers
  !check primary root, with and without vertical distribution
  implicit none
  integer , intent(in) :: I,J,NZ
  character(len=*), intent(in) :: info
  integer :: NR,L,N,NE
  real(r8) :: dRootMyco1stElm_raxs(NumPlantChemElms)

  associate(                                                         &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,          &
    MY                        => plt_morph%MY,                       &
    NU                        => plt_site%NU,                        &
    NL                        => plt_site%NL,                        &    
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr, &
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,       &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft          &
  )
  DO NR=1,NumRootAxes_pft(NZ)
    DO N=1,MY(NZ)
      dRootMyco1stElm_raxs(:)=0._r8
      DO L=NU,NL
        DO NE=1,NumPlantChemElms
          dRootMyco1stElm_raxs(NE)=dRootMyco1stElm_raxs(NE)+RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)  
        ENDDO
      ENDDO
      write(111,*)trim(info),' rootcheck:',I+J/24.,NZ,N,NR,&
         dRootMyco1stElm_raxs(:),RootMyco1stElm_raxs(:,N,NR,NZ),dRootMyco1stElm_raxs(:)-RootMyco1stElm_raxs(:,N,NR,NZ)
      dRootMyco1stElm_raxs(:)=dRootMyco1stElm_raxs(:)-RootMyco1stElm_raxs(:,N,NR,NZ)   
      if(any(abs(dRootMyco1stElm_raxs(:))>1.e-8_r8))stop
    ENDDO
  ENDDO
  end associate
  end subroutine RootCheck
end module RootMod
