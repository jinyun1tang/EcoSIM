module RootMod
  use data_kind_mod,       only: r8 => DAT_KIND_R8,yearIJ_type
  use minimathmod,         only: safe_adb, AZMAX1, AZMIN1, AZERO
  use EcoSIMCtrlMod,       only: lsoilCompaction
  use DebugToolMod,        only: PrintInfo
  use PlantNonstElmDynMod, only: RepleteLowSeaStorByRoot
  use PlantBalMod,         only: SumRootBiome
  use EcosimConst
  use PlantBGCPars
  use ElmIDMod
  use PlantMathFuncMod
  use PlantAPIData
  use NoduleBGCMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: RootBGCModel
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine RootBGCModel(yearIJ,NZ,TFN6_vr,CNRTW,CPRTW,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)
!
! Roots and mycorrhizae both have primary and secondary components
! primary roots are accordant to shoots, seminar root corresponds to main branch.
! not all plants have seminal roots.
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer , intent(in) :: NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)                        !temperature function for root maintenance respiration
  real(r8), intent(in) :: CNRTW,CPRTW                         !NC, PC mass ratio of root growth
  real(r8), intent(out) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(out) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(out) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(out) :: RootSinkC(pltpar%jroots)

  character(len=*), parameter :: subname='RootBGCModel'
  integer, parameter  :: NumRootAxes4DeadPlant =0    !
  real(r8) :: TotAbsorbRootVol
  real(r8) :: litrflx(NumPlantChemElms)
  logical :: flag2ndGrowth(JZ1,MaxNumRootAxes)  
  real(r8) :: RCO2flx,err,resp0
  integer :: N,M,NE,L,K
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)

  associate(                                                            &
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft       ,& !input  :plant stored nonstructural element at current step, [g d-2]
    ZERO4LeafVar_pft           => plt_biom%ZERO4LeafVar_pft            ,& !input  :threshold zero for leaf calculation, [-]
    DLYR3                      => plt_site%DLYR3                       ,& !input  :vertical thickness of soil layer, [m]
    PlantPopulation_pft        => plt_site%PlantPopulation_pft         ,& !input  :plant population, [d-2]
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    NU                         => plt_site%NU                          ,& !input  :current soil surface layer number, [-]    
    MaxSoiL4Root_pft           => plt_morph%MaxSoiL4Root_pft           ,& !input  :maximum soil layer number for all root axes,[-]        
    Myco_pft                   => plt_morph%Myco_pft                   ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolPattern_pft    => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr            ,& !input  :plant LitrFall element, [g d-2 h-1]        
    RootNutUptakeN_pft         => plt_rbgc%RootNutUptakeN_pft          ,& !inoput :total N uptake by plant roots, [gN d-h2 h-1]
    RootNutUptakeP_pft         => plt_rbgc%RootNutUptakeP_pft          ,& !inoput :total P uptake by plant roots, [gP d-h2 h-1]
    SeedMeanLen_pft            => plt_morph%SeedMeanLen_pft            ,& !input  :seed length, [m]
    RootVH2O_pvr               => plt_morph%RootVH2O_pvr               ,& !input  :root layer volume water, [m2 d-2]
    RootCO2Autor_pvr           => plt_rbgc%RootCO2Autor_pvr            ,& !input  :root respiration constrained by O2, [g d-2 h-1]   
    NMaxRootBotLayer_pft       => plt_morph%NMaxRootBotLayer_pft       ,& !input  :maximum soil layer number for all root axes, [-]         
    RootSAreaPerPlant_pvr      => plt_morph%RootSAreaPerPlant_pvr      ,& !input  :root layer area per plant, [m p-1]
    RootPoreVol_pvr            => plt_morph%RootPoreVol_pvr            ,& !input  :root layer volume air, [m2 d-2]
    RootLenDensPerPlant_pvr    => plt_morph%RootLenDensPerPlant_pvr    ,& !input  :root layer length density, [m m-3]
    RootTotLenPerPlant_pvr     => plt_morph%RootTotLenPerPlant_pvr     ,& !inoput :root layer length per plant, [m p-1]
    RootAbsorbLenPerPlant_pvr  => plt_morph%RootAbsorbLenPerPlant_pvr  ,& !inoput :root layer length per plant, [m p-1]    
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft         ,& !input  :soil layer at planting depth, [-]
    RootPorosity_pft           => plt_morph%RootPorosity_pft           ,& !input  :root porosity, [m3 m-3]
    SeedVolumeMean_pft         => plt_morph%SeedVolumeMean_pft         ,& !input  :seed volume, [m3 ]
    SeedAreaMean_pft           => plt_morph%SeedAreaMean_pft           ,& !input  :seed surface area, [m2]
    NumPrimeRootAxes_pft       => plt_morph%NumPrimeRootAxes_pft       ,& !input  :root primary axis number,[-]
    iPlantRootState_pft        => plt_pheno%iPlantRootState_pft        ,& !output :flag to detect root system death,[-]
    iPlantShootState_pft       => plt_pheno%iPlantShootState_pft        & !output :flag to detect canopy death,[-]
  )
!     ROOT GROWTH
!
  call PrintInfo('beg '//subname)
  !  call RootCheck(I,J,NZ,'head')

  !first add newly acquired nutrients
  call SummarizeRootSink(yearIJ%I,yearIJ%J,NZ,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,flag2ndGrowth)

!  resp0=sum(plt_rbgc%RootCO2Autor_pvr(1:Myco_pft(NZ),NU:NMaxRootBotLayer_pft(NZ),NZ))
!  write(444,*)yearIJ%I*1000+yearIJ%J/24.,'beg',resp0,'maxl=',NMaxRootBotLayer_pft(NZ)

!  call SumRootBiome(yearIJ,NZ,mass_inital)

  call RootBiochemistry(yearIJ%I,yearIJ%J,NZ,TFN6_vr,CNRTW,CPRTW,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,flag2ndGrowth)
 
  !  call RootCheck(I,J,NZ,'rootbiochem')
    
  !
  !     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
  !     GERMINATION)
  !
!  call SumRootBiome(yearIJ,NZ,mass_finale)
!  write(444,*)yearIJ%I*1000+yearIJ%J/24.,'end',mass_inital(ielmc),mass_finale(ielmc),mass_inital(ielmc)-mass_finale(ielmc), &
!    sum(plt_rbgc%RootCO2Autor_pvr(1:Myco_pft(NZ),NU:NMaxRootBotLayer_pft(NZ),NZ))-resp0
  
  RootTotLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)    = RootTotLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+SeedMeanLen_pft(NZ)
  RootAbsorbLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootAbsorbLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+SeedMeanLen_pft(NZ)
  
  IF(DLYR3(NGTopRootLayer_pft(NZ)).GT.ZERO)THEN
    !by m3
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)= &
      RootTotLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)/DLYR3(NGTopRootLayer_pft(NZ))
  ELSE
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=0._r8
  ENDIF

  TotAbsorbRootVol=RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedVolumeMean_pft(NZ)*PlantPopulation_pft(NZ)
  RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)       = RootPorosity_pft(ipltroot,NZ)*TotAbsorbRootVol
  RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)          = (1.0_r8-RootPorosity_pft(ipltroot,NZ))*TotAbsorbRootVol
  RootSAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootSAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+&
    SeedAreaMean_pft(NZ)

  IF(NumRootAxes4DeadPlant.EQ.NumPrimeRootAxes_pft(NZ) .OR. (SeasonalNonstElms_pft(ielmc,NZ).LE.ZERO4LeafVar_pft(NZ).AND. &
    iPlantPhenolPattern_pft(NZ).NE.iplt_annual))THEN
    iPlantRootState_pft(NZ)  = iDead
    iPlantShootState_pft(NZ) = iDead
  ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
  call RootNodulBiochemistry(yearIJ%I,yearIJ%J,NZ,TFN6_vr)

!   call SumRootBiome(NZ,mass_finale)
  call PrintInfo('end '//subname)  
  end associate
  end subroutine RootBGCModel

!----------------------------------------------------------------------------------------------------
  subroutine RootBiochemistry(I,J,NZ,TFN6_vr,CNRTW,CPRTW,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,flag2ndGrowth)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: CNRTW,CPRTW    !NC, PC mass ratio of root growth
  REAL(R8), INTENT(in) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), INTENT(in) :: RootSinkC(pltpar%jroots)
  real(r8), INTENT(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  logical, intent(in) :: flag2ndGrowth(JZ1,MaxNumRootAxes)
  real(r8) :: TotPopuRoot2ndLlenAxes
  integer :: LL,LZ,Lnext,L,N
  real(r8) :: TotRoot2ndPopuC !primary/secondary root C
  logical  :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)  !Primary root tip zone update flag
  logical  :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: masst_inital(NumPlantChemElms)
  real(r8) :: masst_finale(NumPlantChemElms)  
  real(r8) :: litrflxt(NumPlantChemElms)
  real(r8) :: RCO2flxt
  character(len=*), parameter :: subname='RootBiochemistry'

!     begin_execution
  associate(                                                                 &
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft         ,& !input  :plant growth type (vascular, non-vascular),[-]
    PSIRoot_pvr                  => plt_ew%PSIRoot_pvr                      ,& !input  :root total water potential, [Mpa]
    PSIRootTurg_vr               => plt_ew%PSIRootTurg_vr                   ,& !input  :root turgor water potential, [Mpa]
    VLSoilPoreMicP_vr            => plt_soilchem%VLSoilPoreMicP_vr          ,& !input  :volume of soil layer, [m3 d-2]
    NU                           => plt_site%NU                             ,& !input  :current soil surface layer number, [-]
    MaxNumRootLays               => plt_site%MaxNumRootLays                 ,& !input  :maximum root layer number,[-]
    ZEROS2                       => plt_site%ZEROS2                         ,& !input  :threshold zero for numerical stability,[-]
    NL                           => plt_site%NL                             ,& !input  :lowest soil layer number,[-]
    Myco_pft                     => plt_morph%Myco_pft                      ,& !input  :mycorrhizal type (no or yes),[-]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft            ,& !input  :soil layer at planting depth, [-]
    fRootGrowPSISense_pvr        => plt_pheno%fRootGrowPSISense_pvr         ,& !output :water stress to plant root growth, [-]
    NMaxRootBotLayer_pft         => plt_morph%NMaxRootBotLayer_pft           & !output :maximum soil layer number for all root axes, [-]
  )
  call PrintInfo('beg '//subname)
  Root1stTipUpdateFlag     = .false. !initialize all primary roots as not updated
  FoundRootAxesTip        = .false.
  NMaxRootBotLayer_pft(NZ) = NGTopRootLayer_pft(NZ) !initialize diagnostic for the deepest root layer
  !
  !     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
  !
  litrflxt=0._r8;RCO2flxt=0._r8
  CALL GetPlantRoot1stDepz(I,J,NZ)

!  call FindMatureRootSegs(I,J,NZ)

  D5010: DO N=1,Myco_pft(NZ)
    D5000: DO L=NU,MaxNumRootLays
      !
      !     IDENTIFY NEXT LOWER ROOT LAYER
      !
      !     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
      !
      IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
        D5003: DO LZ=L+1,NL
          IF(VLSoilPoreMicP_vr(LZ).GT.ZEROS2 .OR. LZ.EQ.NL)THEN
            Lnext=LZ
            EXIT
          ENDIF
        ENDDO D5003

      !   call SumRootBiome(NZ,masst_inital)

        IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
          fRootGrowPSISense_pvr(N,L,NZ)=EXP(0.05_r8*AMAX1(PSIRoot_pvr(N,L,NZ),-5000._r8))
        ELSE
          fRootGrowPSISense_pvr(N,L,NZ)=EXP(0.10_r8*AMAX1(PSIRoot_pvr(N,L,NZ),-5000._r8))
        ENDIF  
        
        call GrowRootMycoAxes(I,J,N,L,Lnext,NZ,FoundRootAxesTip,Root1stTipUpdateFlag,TFN6_vr,&
          RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,flag2ndGrowth,CNRTW,CPRTW,&
          fRootGrowPSISense_pvr(N,L,NZ),TotPopuRoot2ndLlenAxes,TotRoot2ndPopuC,litrflxt,RCO2flxt)
        
      !   call SumRootBiome(NZ,masst_finale)

!====================================================================================
        !masst_inital=masst_finale
        !     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
        !     SEASONAL STORAGE POOL IS DEPLETED
        call RepleteLowSeaStorByRoot(I,J,N,L,NZ)

        call DiagRootGeometry(I,J,N,L,NZ,flag2ndGrowth,TotRoot2ndPopuC,TotPopuRoot2ndLlenAxes)
      !   call SumRootBiome(NZ,masst_finale)
      ENDIF
    ENDDO D5000
  ENDDO D5010

  call PrintInfo('end '//subname)
  end associate
  end subroutine RootBiochemistry
!----------------------------------------------------------------------------------------------------

  subroutine DiagRootGeometry(I,J,N,L,NZ,flag2ndGrowth,TotRoot2ndPopuC,TotPopuRoot2ndLlenAxes)
  implicit none
  integer , intent(in)  :: I,J,N,L,NZ
  logical,  intent(in)  :: flag2ndGrowth(JZ1,MaxNumRootAxes)    !flag for primary roots that are undergoing secondary growth  
  real(r8), intent(in)  :: TotRoot2ndPopuC                      !total secondary root C mass
  real(r8), intent(inout)  :: TotPopuRoot2ndLlenAxes            !total secondary root length  
  real(r8) :: CoarseVol                    !coarse root volume, [m3]    
  real(r8) :: TotRoot1stPopuC              !total primary root C mass
  real(r8) :: TotRoot1stLlenAxesPP         !total primary root length  
  real(r8) :: TotPopuRoot1stLen 
  real(r8) :: TotPopuRootLen               !total root length
  real(r8) :: TotFineRootC                 !total root C mass
  real(r8) :: TotAbsorbRootVol             !total volume
  real(r8) :: Root1stSurfAbsorbArea,Root2ndSurfArea
  real(r8) :: Root1stRadiusEst,SCAL,FineRoot1stTotc,TotFineRoot1stLlenAxesPP  
  real(r8) :: TotCoarseRoot1stLlenAxesPP,CoarseRoot1stTotC
  real(r8) :: AreaTranspt,DTransptTube,TotAbsorbRootVolPP       
  integer  :: NTG,NE,NR
  logical :: checkRootExt,check2ndGrothRoot
  character(len=*), parameter :: subname='DiagRootGeometry'
  associate(                                                              &  
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft                 ,& !input  :threshold zero for plang growth calculation, [-]  
    PlantPopulation_pft       => plt_site%PlantPopulation_pft            ,& !input  :plant population, [d-2]
    ZERO                      => plt_site%ZERO                           ,& !input  :threshold zero for numerical stability, [-]
    RootPorosity_pft          => plt_morph%RootPorosity_pft              ,& !input  :root porosity, [m3 m-3]
    DLYR3                     => plt_site%DLYR3                          ,& !input  :vertical thickness of soil layer, [m]
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr                      ,& !input  :root total water potential, [Mpa]
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr                   ,& !input  :root turgor water potential, [Mpa]
    Root1stMaxRadius1_pft     => plt_morph%Root1stMaxRadius1_pft         ,& !input  :root radius primary axes, [m]
    Root2ndXSecArea_pft       => plt_morph%Root2ndXSecArea_pft           ,& !input  :root cross-sectional area of secondary axes, [m2]    
    Root2ndMaxRadius1_pft     => plt_morph%Root2ndMaxRadius1_pft         ,& !input  :root radius secondary axes, [m]
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft          ,& !input  :maximum radius of secondary roots, [m]
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft          ,& !input  :maximum radius of primary roots, [m]
    k_fine_comp               => pltpar%k_fine_comp                      ,& !input  :fine litter complex id    
    Root1stXSecArea_pft       => plt_morph%Root1stXSecArea_pft           ,& !input  :root cross-sectional area primary axes, [m2]  
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft         ,& !input  :plant growth type (vascular, non-vascular),[-]    
    RootMatureAge_pft         => plt_morph%RootMatureAge_pft             ,& !input : Root maturation age, [h]
    NumPrimeRootAxes_pft      => plt_morph%NumPrimeRootAxes_pft          ,& !input  :root primary axis number,[-]    
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm             ,& !input  :C woody/soft fraction in root,[-]    
    Root1stLenPP_rpvr         => plt_morph%Root1stLenPP_rpvr             ,& !input  :root layer length primary axes, [m d-2]    
    Root2ndXNumL_rpvr         => plt_morph%Root2ndXNumL_rpvr             ,& !input  :root layer number axes, [d-2]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr      ,& !input  :root layer element primary axes, [g d-2]    
    RootAge_rpvr              => plt_morph%RootAge_rpvr                  ,& !inoput :root age,[h]
    CoarseRootVolPerMassC_pft => plt_morph%CoarseRootVolPerMassC_pft     ,& !input  :coarse root volume:mass ratio, [m3 gC-1]
    FineRootVolPerMassC_pft   => plt_morph%FineRootVolPerMassC_pft       ,& !input  :fine root volume:mass ratio, [m3 gC-1]
    Root1stRadius_rpvr        => plt_morph%Root1stRadius_rpvr            ,& !input  :root layer radius for each primary axes, [m]    
    Root1stTransptArea_pvr    => plt_morph%Root1stTransptArea_pvr        ,& !output :mean transport area by 1st order root, [m2 d-2]         
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr             ,& !output :root layer radius primary axes, [m]    
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr                  ,& !output :root space volume occupied by water in each layer, [m2 d-2]   
    Root2ndRadius_rpvr        => plt_morph%Root2ndRadius_rpvr            ,& !output :root layer radius secondary axes, [m]     
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr       ,& !output :layer root length density, [m m-3]
    RootArea1stPP_pvr         => plt_morph%RootArea1stPP_pvr             ,& !output :layer 1st root area per plant, [m2 plant-1]
    RootArea2ndPP_pvr         => plt_morph%RootArea2ndPP_pvr             ,& !output :layer 2nd root area per plant, [m2 plant-1] 
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr               ,& !output :layer root volume air, [m2 d-2]
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr                ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr                ,& !inoput :root aqueous content, [g d-2]
    RootSAreaPerPlant_pvr     => plt_morph%RootSAreaPerPlant_pvr         ,& !output :layer 1st+2nd root area per plant, [m2 plant-1]    
    RootTotLenPerPlant_pvr    => plt_morph%RootTotLenPerPlant_pvr        ,& !output :total root length per plant in layer (including root hair), [m p-1]
    RootAbsorbLenPerPlant_pvr => plt_morph%RootAbsorbLenPerPlant_pvr     ,& !output :total absorptive root length per plant in layer, [m p-1]
    Root2ndEffLen4uptk_rpvr   => plt_morph%Root2ndEffLen4uptk_rpvr       ,& !output :root layer average length, [m]
    RootLenPerPlant_pvr       => plt_morph%RootLenPerPlant_pvr           ,& !output :root length per plant in layer (excluding root hair), [m p-1]    
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft          & !inoput :gaseous flux fron root disturbance, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)
  !
  CoarseRoot1stTotC              = 0._r8
  FineRoot1stTotc                = 0._r8
  TotFineRoot1stLlenAxesPP       = 0._r8
  TotCoarseRoot1stLlenAxesPP     = 0._r8
  TotAbsorbRootVolPP             = 0._r8
  Root1stSurfAbsorbArea          = 0._r8
  Root1stTransptArea_pvr(N,L,NZ) = 0._r8

  DO NR=1,NumPrimeRootAxes_pft(NZ)
    ! TotRoot1stLlenAxesPP=total primary root length per plant
    ! TotRoot1stPopuC=total primary root C mass for whole population     
    if(N==ipltroot)then
      DTransptTube     = AMIN1(ZSTX,AMAX1(FSTK*Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(N,NZ)))
      Root1stRadiusEst = AMAX1(Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
      !transport area for axis NR, 
      AreaTranspt      = PICON*(2._r8*Root1stRadiusEst*DTransptTube-DTransptTube**2)

      if(flag2ndGrowth(L,NR) .and. Root1stLenPP_rpvr(L,NR,NZ)>0._r8)then
        CoarseRoot1stTotC          = CoarseRoot1stTotC+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)
        TotCoarseRoot1stLlenAxesPP = TotCoarseRoot1stLlenAxesPP+Root1stLenPP_rpvr(L,NR,NZ)
      else
        FineRoot1stTotc          = FineRoot1stTotc+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)
        TotFineRoot1stLlenAxesPP = TotFineRoot1stLlenAxesPP+Root1stLenPP_rpvr(L,NR,NZ)
        Root1stSurfAbsorbArea    = Root1stSurfAbsorbArea+TwoPiCON*Root1stRadiusEst*Root1stLenPP_rpvr(L,NR,NZ)
        TotAbsorbRootVolPP       = TotAbsorbRootVolPP+AreaTranspt*Root1stLenPP_rpvr(L,NR,NZ)
      endif      
      Root1stTransptArea_pvr(N,L,NZ)=Root1stTransptArea_pvr(N,L,NZ)+AreaTranspt
    endif              
  ENDDO
  if(NumPrimeRootAxes_pft(NZ)>0)Root1stTransptArea_pvr(N,L,NZ)=Root1stTransptArea_pvr(N,L,NZ)/NumPrimeRootAxes_pft(NZ)

  ! ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
  ! TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
  ! account for fine roots only

  IF(N.EQ.ipltroot)THEN
    TotRoot1stLlenAxesPP   = TotCoarseRoot1stLlenAxesPP+TotFineRoot1stLlenAxesPP
    CoarseVol              = CoarseRoot1stTotC*CoarseRootVolPerMassC_pft(NZ)
  else
    TotRoot1stLlenAxesPP = TotCoarseRoot1stLlenAxesPP+TotFineRoot1stLlenAxesPP    
    CoarseVol=0._r8
  ENDIF

  TotFineRootC         = FineRoot1stTotc+TotRoot2ndPopuC
  TotPopuRoot1stLen    = TotRoot1stLlenAxesPP*PlantPopulation_pft(NZ)
  TotPopuRootLen       = TotPopuRoot1stLen+TotPopuRoot2ndLlenAxes

  checkRootExt      = TotPopuRootLen.GT.ZERO4Groth_pft(NZ) .AND. TotFineRootC.GT.ZERO4Groth_pft(NZ) &
    .AND. PlantPopulation_pft(NZ).GT.ZERO4Groth_pft(NZ)

  IF(checkRootExt)THEN  
    !there are roots
    RootAbsorbLenPerPlant_pvr(N,L,NZ) = TotFineRoot1stLlenAxesPP+TotPopuRoot2ndLlenAxes/PlantPopulation_pft(NZ)
    RootTotLenPerPlant_pvr(N,L,NZ) = TotPopuRootLen/PlantPopulation_pft(NZ)
    RootLenPerPlant_pvr(N,L,NZ)    = TotRoot1stLlenAxesPP        
    IF(DLYR3(L).GT.ZERO)THEN
      !per volume
      RootLenDensPerPlant_pvr(N,L,NZ) = RootTotLenPerPlant_pvr(N,L,NZ)/DLYR3(L)
    ELSE
      RootLenDensPerPlant_pvr(N,L,NZ) = 0._r8
    ENDIF
    !
    !total root volume [m3], with FineRootVolPerMassC_pft including contribution from water, fine roots volume > 80% as water
    !coarse roots, 45-60% volume as water 

    check2ndGrothRoot=N.EQ.ipltroot .and. is_plant_treelike(iPlantRootProfile_pft(NZ)) &
      .and. TotPopuRoot1stLen>0._r8 .and. CoarseVol > 0._r8

    
    if(check2ndGrothRoot)then

      !assuming elongation zone is negligible in length, pi*r^2*l=vol/pop       
      Root1stRadiusEst = sqrt(CoarseVol/(PICON*TotPopuRoot1stLen))
      !total absorption root volume
      TotAbsorbRootVol = TotAbsorbRootVolPP*PlantPopulation_pft(NZ) + AMAX1(Root2ndXSecArea_pft(N,NZ)*TotPopuRoot2ndLlenAxes, &
        TotRoot2ndPopuC*FineRootVolPerMassC_pft(N,NZ)*PSIRootTurg_vr(N,L,NZ))
      
      Root1stRadius_pvr(N,L,NZ) = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ),Root1stRadiusEst)
      
    else
      TotAbsorbRootVol=AMAX1(Root1stXSecArea_pft(N,NZ)*TotPopuRoot1stLen+Root2ndXSecArea_pft(N,NZ)*TotPopuRoot2ndLlenAxes &
        ,TotFineRootC*FineRootVolPerMassC_pft(N,NZ)*PSIRootTurg_vr(N,L,NZ))    
      Root1stRadius_pvr(N,L,NZ) = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
      Root1stSurfAbsorbArea     = TwoPiCON*Root1stRadius_pvr(N,L,NZ)*TotPopuRoot1stLen
    endif  

    RootPoreVol_pvr(N,L,NZ) = RootPorosity_pft(N,NZ)*TotAbsorbRootVol
    RootVH2O_pvr(N,L,NZ)    = (1.0_r8-RootPorosity_pft(N,NZ))*TotAbsorbRootVol

    !primary and secondary root radius
    Root2ndRadius_rpvr(N,L,NZ) = AMAX1(Root2ndMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root2ndMaxRadius_pft(N,NZ))
    Root2ndSurfArea            = TwoPiCON*Root2ndRadius_rpvr(N,L,NZ)*TotPopuRoot2ndLlenAxes
    RootArea1stPP_pvr(N,L,NZ)  = Root1stSurfAbsorbArea/PlantPopulation_pft(NZ)
    RootArea2ndPP_pvr(N,L,NZ)  = Root2ndSurfArea/PlantPopulation_pft(NZ)

    IF(Root2ndXNumL_rpvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      Root2ndEffLen4uptk_rpvr(N,L,NZ)=AMAX1(Root2ndTipLen4uptk,TotPopuRoot2ndLlenAxes/Root2ndXNumL_rpvr(N,L,NZ))
    ELSE
      Root2ndEffLen4uptk_rpvr(N,L,NZ)=Root2ndTipLen4uptk
    ENDIF
    RootSAreaPerPlant_pvr(N,L,NZ)=RootArea1stPP_pvr(N,L,NZ)+RootArea2ndPP_pvr(N,L,NZ)
  ELSE
    !no roots
    RootArea1stPP_pvr(N,L,NZ)       = 0._r8
    RootArea2ndPP_pvr(N,L,NZ)       = 0._r8
    RootLenPerPlant_pvr(N,L,NZ)     = 0._r8
    RootTotLenPerPlant_pvr(N,L,NZ)  = 0._r8
    RootLenDensPerPlant_pvr(N,L,NZ) = 0._r8
    RootPoreVol_pvr(N,L,NZ)        = 0._r8
    RootVH2O_pvr(N,L,NZ)            = 0._r8
    RootSAreaPerPlant_pvr(N,L,NZ)    = 0._r8    
    Root1stRadius_pvr(N,L,NZ)       = Root1stMaxRadius_pft(N,NZ)
    Root2ndRadius_rpvr(N,L,NZ)      = Root2ndMaxRadius_pft(N,NZ)
    Root1stTransptArea_pvr(N,L,NZ)  = PICON*Root1stMaxRadius_pft(N,NZ)**2
    Root2ndEffLen4uptk_rpvr(N,L,NZ) = Root2ndTipLen4uptk
    DO NTG=idg_beg,idg_NH3
      RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-(trcg_rootml_pvr(NTG,N,L,NZ)+trcs_rootml_pvr(NTG,N,L,NZ))
    ENDDO
    trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
    trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
  ENDIF

  call PrintInfo('end '//subname)      
  end associate      
  end subroutine DiagRootGeometry
!----------------------------------------------------------------------------------------------------
  subroutine Grow2ndRootAxes(I,J,L,NZ,N,NR,Nutstress4GrossResp,CNRTW,CPRTW,dLext2nd,&
    RespElongWatSens,DMRespEff,fRootGrowPSISense,TFN6_vr,RootSinkC_vr,Root2ndSink_pvr,&
    RootNsink2nd,litrflx2,RCO2flx2,Root2ndStrutRemob)
  !
  !Description:
  !Grow secondary roots
    
  implicit none
  integer, intent(in) :: I,J,L,NZ
  integer, intent(in) :: N                !root type id
  integer, intent(in) :: NR               !root axis id
  real(r8), intent(in) :: Nutstress4GrossResp
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: dLext2nd                  !soil resistance allwed root penetration [h m-1]
  real(r8), intent(in) :: RespElongWatSens          !
  real(r8), intent(in) :: DMRespEff                 !respiraiton quotient
  real(r8), intent(in) :: fRootGrowPSISense  
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: RootNsink2nd  
  real(r8), intent(out) :: litrflx2(NumPlantChemElms)      !litterfall flux 
  real(r8), intent(out) :: RCO2flx2                        !CO2 respiraiton flux 
  real(r8), intent(out) :: Root2ndStrutRemob(NumPlantChemElms)      
  character(len=*), parameter :: subname='Grow2ndRootAxes'
  real(r8) :: Rmaint2nd_CO2
  real(r8) :: RNonstCO2_OUltd,RNonstCO2_Oltd  
  REAL(R8) :: RCO2MaintDef2ndStruct_Oltd  !maintenance deficit met up by structural root biomass remobilization under O2 limitation
  real(r8) :: RCO2MaintDef2ndStruct_OUltd !maintenance deficit met up by structural root biomass remobilization without O2 limitation
  real(r8) :: RCO2Nonst4Nassim_Oltd       !respiration for N assimilation ltd by O2
  real(r8) :: RCO2Nonst4Nassim_OUltd      !respiration for N assimilation unltd by O2
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4GrowC_OUltd
  real(r8) :: RootMycoNonst4GrowC_Oltd    
  real(r8) :: Root2ndPopExtenz  !whole plant population 2ndary root extension
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms)  !nonstrucal biomass converted into biomass 
  real(r8) :: Root2ndNetGrowthElms(NumPlantChemElms)
  real(r8) :: dmass(NumPlantChemElms)  
  real(r8) :: RCO2XMaint2nd_Oltd
  real(r8) :: RGrowCO2_Oltd,RCO2T2nd_Oltd
  real(r8) :: RCO2XMaint2nd_OUltd
  real(r8) :: RGrowCO2_OUltd
  real(r8) :: RCO2T2nd_OUltd
  real(r8) :: DMRTR
  real(r8) :: PPOOLB,ZPOOLB
  real(r8) :: RTN2X,RTN2Y  
  real(r8) :: respscal,Root2ndExtPot  
  real(r8) :: FNP,ratio
  real(r8) :: FracRoot2ndCSinkL,FracRoot2ndSinkN,FracRoot2ndSinkP  
  integer  :: M,NE,NumGroFineAxes

  associate(                                                          &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr     ,& !input  :O2 constraint to root respiration (0-1), [-]
    rProteinC2RootN_pft       => plt_allom%rProteinC2RootN_pft       ,& !input  :Protein C to leaf N ratio in roots, [-]
    rProteinC2RootP_pft       => plt_allom%rProteinC2RootP_pft       ,& !input  :Protein C to leaf P ratio in roots, [-]
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm         ,& !input  :C woody/soft fraction in root,[-]
    CNRTS_pft                 => plt_allom%CNRTS_pft                 ,& !input  :root N:C ratio x root growth yield, [-]
    CPRTS_pft                 => plt_allom%CPRTS_pft                 ,& !input  :root P:C ratio x root growth yield, [-]
    RootBiomGrosYld_pft       => plt_allom%RootBiomGrosYld_pft       ,& !input  :root growth yield, [g g-1]
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft     ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft      ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    fTgrowRootP_vr            => plt_pheno%fTgrowRootP_vr            ,& !input  :root layer temperature growth functiom, [-]
    DLYR3                     => plt_site%DLYR3                      ,& !input  :vertical thickness of soil layer, [m]
    RootBranchFreq_pft        => plt_morph%RootBranchFreq_pft        ,& !input  :root brancing frequency, [m-1]
    Root2ndSpecLen_pft        => plt_morph%Root2ndSpecLen_pft        ,& !input  :specific root length secondary axes, [m gC-1]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft         ,& !input  :number of main branch,[-]
    GrainFillDowreg_brch      => plt_photo%GrainFillDowreg_brch      ,& !input  :grain fill down-regulation of annual plants, [-]
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !inoput :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr           ,& !inoput :root layer protein C, [gC d-2]
    RootCO2Autor_pvr          => plt_rbgc%RootCO2Autor_pvr           ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr         ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr         ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootMaintDef_CO2_pvr      => plt_bgcr%RootMaintDef_CO2_pvr       ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]
    Root2ndXNumL_rpvr         => plt_morph%Root2ndXNumL_rpvr         ,& !inoput :root layer number axes, [d-2]
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr           ,& !inoput :root layer length secondary axes, [m d-2]
    NumAxesPerPrimRoot_pft    => plt_morph%NumAxesPerPrimRoot_pft    ,& !output :primary root axes number, [d-2]            
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr           & !output :root layer number secondary axes, [d-2]
  )
  !     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
  !     TO CURRENT AXIS
  !
  !     Root2ndSink_pvr=total secondary root sink strength
  !     RootSinkC_vr=total root sink strength
  !     FracRoot2ndCSinkL=fraction of secondary root sink strength in axis
  !
  call PrintInfo('beg '//subname)
  DO NE=1,NumPlantChemElms
    dmass(NE)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)
  ENDDO
  
  litrflx2=0._r8
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot2ndCSinkL=Root2ndSink_pvr(N,L,NR)/RootSinkC_vr(N,L)
  ELSE
    FracRoot2ndCSinkL=1.0_r8
  ENDIF
  
  !
  !     N,P CONSTRAINT ON SECONDARY ROOT RESPIRATION FROM
  !     NON-STRUCTURAL C:N:P
  !
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     Nutstress4GrossResp=N,P constraint on growth respiration
  !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !

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

  !if herbaceous root or drought deciduous, maintenance is moisture dependent.
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)) .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint2nd_CO2=Rmaint2nd_CO2*fRootGrowPSISense
  ENDIF
!
!     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RNonstCO2_OUltd=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     
!     FracRoot2ndCSinkL=fraction of secondary root sink strength in axis
!     CPOOL=non-structural C mass
!     fTgrowRootP_vr=temperature function for root growth
!     Nutstress4GrossResp=N,P constraint on respiration
!     GrainFillDowreg_brch=termination feedback inhibition on C3 CO2
!     fRootGrowPSISense=growth function of root water potential
!
  respscal=VMXC*FracRoot2ndCSinkL*fTgrowRootP_vr(L,NZ)*fRootGrowPSISense &
    *Nutstress4GrossResp*GrainFillDowreg_brch(MainBranchNum_pft(NZ),NZ)

  RNonstCO2_OUltd=AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))*respscal

!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RNonstCO2_Oltd=respiration from non-structural C limited by O2
!     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!     RCO2XMaint2nd_OUltd,RCO2XMaint2nd_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
!     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration unltd,ltd by O2 and unlimited by N,P
!     RespElongWatSens=respiration function of root water potential
!
  RNonstCO2_Oltd               = RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RCO2XMaint2nd_OUltd          = RNonstCO2_OUltd-Rmaint2nd_CO2
  RCO2XMaint2nd_Oltd           = RNonstCO2_Oltd-Rmaint2nd_CO2
  RGrowCO2_OUltd               = AZMAX1(RCO2XMaint2nd_OUltd)*RespElongWatSens
  RGrowCO2_Oltd                = AZMAX1(RCO2XMaint2nd_Oltd)*RespElongWatSens
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RCO2XMaint2nd_Oltd,0._r8)
  !
  !     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
  !     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
  !
  !     FracRoot2ndCSinkL=fraction of secondary root sink strength in axis
  !     FNP=growth respiration limited by non-structural N,P
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !
  DMRTR  = DMRespEff
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*FracRoot2ndCSinkL
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*FracRoot2ndCSinkL
  FNP    = AMIN1(ZPOOLB/CNRTS_pft(NZ),PPOOLB/CPRTS_pft(NZ))*DMRTR

  IF(RGrowCO2_OUltd.GT.0.0_r8)THEN
    RGrowCO2_OUltd=AMIN1(RGrowCO2_OUltd,FNP)
  ELSE
    RGrowCO2_OUltd=0._r8
  ENDIF

  !O2-limited active growth
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
!     DMRespEff=root C respiration vs nonstructural C consumption, i.e. respiration quotient
!     RootBiomGrosYld_pft=root growth yield
!
  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRespEff
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRespEff

  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)
  RCO2Nonst4Nassim_OUltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_OUltd(ielmn))

  !O2-limited biomass growth
  RootMycoNonst4Grow_Oltd(ielmc) = RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn) = AZMAX1(AMIN1(FracRoot2ndCSinkL*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
  RootMycoNonst4Grow_Oltd(ielmp) = AZMAX1(AMIN1(FracRoot2ndCSinkL*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))  
  RCO2Nonst4Nassim_Oltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_Oltd(ielmn))


  ! apply turgor pressure limitation on elongation
  Root2ndExtPot    = RootMycoNonst4Grow_Oltd(ielmc)*FracRootElmAllocm(ielmc,k_fine_comp)*Root2ndSpecLen_pft(N,NZ)
  Root2ndPopExtenz = Root2ndExtPot

  if(lsoilCompaction)then
    if(Root2ndExtPot>0._r8)then  
      Root2ndPopExtenz = AMIN1(Root2ndExtPot,dLext2nd)
      if(Root2ndPopExtenz <Root2ndExtPot)then    
        ratio=Root2ndPopExtenz/Root2ndExtPot
        DO NE=1,NumPlantChemElms
          RootMycoNonst4Grow_Oltd(NE)=RootMycoNonst4Grow_Oltd(NE)*ratio
        ENDDO
      endif
    else
      Root2ndPopExtenz=0._r8
    endif
  endif

  !Consume nonstrucal elements for growth
  DO NE=1,NumPlantChemElms
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Grow_Oltd(NE)    
  ENDDO

  ! net growth=positive growth - senescence
  Root2ndNetGrowthElms = RootMycoNonst4Grow_Oltd

  !
  !     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
  !     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
  !     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LitrFall
  !
  call DiagSenes2ndRootAxes(I,J,N,NR,L,NZ,RCO2XMaint2nd_OUltd,RCO2XMaint2nd_Oltd,Root2ndPopExtenz,Root2ndNetGrowthElms,&
    litrflx2,RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd)

  !
  !     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
  !
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
  !     Rmaint2nd_CO2=root maintenance respiration
  !     RNonstCO2_Oltd=respiration from non-structural C limited by O2
  !     RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn ltd by O2
  !     RootMycoNonst4Grow_Oltd(:)=nonstructural N,P ltd by O2 used in growth

  !do final summary
  !
  !     TOTAL SECONDARY ROOT RESPIRATION
  !
  !     RCO2T2nd_OUltd,RCO2T2nd_Oltd=total C respiration unltd,ltd by O2
  !     Rmaint2nd_CO2=root maintenance respiration
  !     RNonstCO2_OUltd,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !     RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd=excess maintenance respiration unltd,ltd by O2
  !     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !     RootCO2Autor_pvr=total root respiration
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  RCO2T2nd_OUltd = AMIN1(Rmaint2nd_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd+RCO2Nonst4Nassim_OUltd+RCO2MaintDef2ndStruct_OUltd
  RCO2T2nd_Oltd  = AMIN1(Rmaint2nd_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+RCO2Nonst4Nassim_Oltd+RCO2MaintDef2ndStruct_Oltd
  RCO2T2nd_Oltd  = AMIN1(RCO2T2nd_Oltd,AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)))
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T2nd_Oltd

  RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T2nd_OUltd
  RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T2nd_Oltd
  RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T2nd_Oltd
  RCO2flx2                   = -RCO2T2nd_Oltd

  !
  !     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
  !
  !     Root2ndPopExtenz=secondary root length extension
  !     RootMycoNonst4Grow_Oltd(ielmc)=secondary root C growth ltd by O2
  !     Root2ndSpecLen_pft=specific secondary root length from startq.f
  !     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
  !     Frac2Senes2=fraction of secondary root C to be remobilized
  !     Root2ndLen_rpvr=secondary root length
  !
  !     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
  !     AND AXIS NUMBER
  !
  !     RootProteinC_pvr=total root protein C mass
  !     RootBranchFreq_pft=root branching frequency from PFT file for generate secondary roots
  !     Root2ndXNum_rpvr,Root2ndXNumL_rpvr=number of secondary root axes
  !
  Root2ndLen_rpvr(N,L,NR,NZ)=Root2ndLen_rpvr(N,L,NR,NZ)+Root2ndPopExtenz
  DO NE=1,NumPlantChemElms
    RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)+Root2ndNetGrowthElms(NE)
  ENDDO
  
  if(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)>0._r8)then  
    RootProteinC_pvr(N,L,NZ)=RootProteinC_pvr(N,L,NZ)+&
      AMIN1(rProteinC2RootN_pft(NZ)*RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ),&
        rProteinC2RootP_pft(NZ)*RootMyco2ndStrutElms_rpvr(ielmp,N,L,NR,NZ))

    !secondary/fine root axes (and root hair) addition is a quadratic function of branching frequency
    RTN2X          = RootBranchFreq_pft(NZ)*NumAxesPerPrimRoot_pft(NZ)    !fine roots
    RTN2Y          = RootBranchFreq_pft(NZ)*RTN2X                       !root hairs
    NumGroFineAxes = (RTN2X+RTN2Y)*DLYR3(L)

    Root2ndXNum_rpvr(N,L,NR,NZ) = NumGroFineAxes
    Root2ndXNumL_rpvr(N,L,NZ)   = Root2ndXNumL_rpvr(N,L,NZ)+Root2ndXNum_rpvr(N,L,NR,NZ)
  endif

  DO NE=1,NumPlantChemElms
    dmass(NE)=dmass(NE)-RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)-litrflx2(NE)
  ENDDO
  dmass(ielmc)=dmass(ielmc)-RCO2T2nd_Oltd
!  write(555,*)I*1000+J/24.,dmass,RCO2T2nd_Oltd
  call PrintInfo('end '//subname)
  end associate  
  end subroutine Grow2ndRootAxes
!----------------------------------------------------------------------------------------------------

  subroutine DiagSenes2ndRootAxes(I,J,N,NR,L,NZ,RCO2XMaint2nd_OUltd,RCO2XMaint2nd_Oltd,Root2ndPopExtenz,&
    Root2ndNetGrowthElms,litrflx2,RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd)
  !
  !Description:
  !Diagnose secondary root senescence and remolization
  implicit none
  integer, intent(in) :: I,J,N,NR,L,NZ
  REAL(R8), intent(in) :: RCO2XMaint2nd_OUltd                        !O2-unlimited growth respiraiton, [gC h-1]
  real(r8), intent(in) :: RCO2XMaint2nd_Oltd                         !O2-limited growth respiraiton, [gC h-1]
  real(r8), intent(out) :: RCO2MaintDef2ndStruct_OUltd               !O2-unlimited maintenance respiraiton deficit [gC h-1]
  real(r8), intent(out) :: RCO2MaintDef2ndStruct_Oltd                !O2-limited maintenance respiraiton deficit [gC h-1]
  real(r8), intent(inout) :: litrflx2(NumPlantChemElms)                !total litterfall flux
  real(r8), intent(inout) :: Root2ndNetGrowthElms(NumPlantChemElms)  !structural biomass growth for 2ndary roots
  real(r8), intent(inout) :: Root2ndPopExtenz                        !potential length extension [m]
  real(r8) :: Root2ndStrutRemob(NumPlantChemElms)                   !remobilization of C,N,P from senescing root
  real(r8) :: Remobl2ndcycl(NumPlantChemElms)                       !recycled nutrients by remobilization
  real(r8) :: Remobl2ndKill(NumPlantChemElms)                       !root lost due to maintenance deficit

  real(r8) :: CCC,CNC,CPC
  real(r8) :: Frac2Senes2         !Fraction of secondary root C to be remobilized
  real(r8) :: RCCE(NumPlantChemElms)
  real(r8) :: RCCC,RCCN,RCCP,FracRemobl  
  real(r8) :: dsenecE,dfineE,dwoodyE
  integer :: M,NE

  associate(                                                          &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]  
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft         ,& !input  :number of main branch,[-]  
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm         ,& !input  :C woody/soft fraction in root,[-]    
    RootNonstructElmConc_rpvr => plt_biom%RootNonstructElmConc_rpvr  ,& !input  :root layer nonstructural C concentration, [g g-1]    
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !input  :root layer element secondary axes, [g d-2]    
    ZERO                      => plt_site%ZERO                       ,& !input  :threshold zero for numerical stability, [-]    
    icwood                    => pltpar%icwood                       ,& !input  :group id of coarse woody litter
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr           ,& !inoput :root layer length secondary axes, [m d-2]    
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]    
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id    
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr     ,& !input  :O2 constraint to root respiration (0-1), [-]    
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]    
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]    
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch        & !input  :plant growth stage, [-]
  )
  !     iPlantCalendar_brch(ipltcal_Emerge),=emergence date
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !     RCCZR,RCCYR=min,max fractions for root C recycling
  !     RCCXR,RCCQR=max fractions for root N,P recycling
  !
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 .AND. &
    RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmn,N,L,NZ) &
      ,RootNonstructElmConc_rpvr(ielmn,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CNKI) &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmp,N,L,NZ) &
      ,RootNonstructElmConc_rpvr(ielmp,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8 &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ) &
      ,RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8 &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ) &
      ,RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+ RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC = 0._r8
    CNC = 0._r8
    CPC = 0._r8
  ENDIF
  RCCC = RCCZR+CCC*RCCYR
  RCCN = CNC*RCCXR
  RCCP = CPC*RCCQR
  RCCE(ielmc)=RCCC
  RCCE(ielmn)=(RCCN+(1.0_r8-RCCN)*RCCC)
  RCCE(ielmp)=(RCCP+(1.0_r8-RCCP)*RCCC)
  !
  !     RECOVERY OF REMOBILIZABLE N,P FROM SECONDARY ROOT DURING
  !     REMOBILIZATION DEPENDS ON ROOT NON-STRUCTURAL C:N:P
  !
  !     RCO2XMaint2nd_OUltd,RCO2XMaint2nd_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
  !     RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd=excess maintenance respiration unltd,ltd by O2
  !     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
  !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !     Root2ndStrutRemob(:)=remobilization of C,N,P from senescing root
  !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !     Frac2Senes2=fraction of secondary root C to be remobilized
  !
  IF(-RCO2XMaint2nd_OUltd.GT.0.0_r8)THEN
    IF(-RCO2XMaint2nd_OUltd.LT.RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
      RCO2MaintDef2ndStruct_OUltd=-RCO2XMaint2nd_OUltd
    ELSE
      RCO2MaintDef2ndStruct_OUltd=AZMAX1(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)
    ENDIF
  ELSE
    RCO2MaintDef2ndStruct_OUltd=0._r8
  ENDIF

  !there is maintenance deficit,
  IF(-RCO2XMaint2nd_Oltd.GT.0.0_r8)THEN
    IF(-RCO2XMaint2nd_Oltd.LT.RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)THEN
      RCO2MaintDef2ndStruct_Oltd=-RCO2XMaint2nd_Oltd
    ELSE
      RCO2MaintDef2ndStruct_Oltd=AZMAX1(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)*RAutoRootO2Limter_rpvr(N,L,NZ)
    ENDIF
  ELSE
    RCO2MaintDef2ndStruct_Oltd=0._r8
  ENDIF

  !Maintenance deficit leads to remobilization/degradation of structural biomass to pay off the deficit
  IF(RCO2MaintDef2ndStruct_Oltd.GT.0.0_r8 .AND. RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    DO NE=1,NumPlantChemElms
      Root2ndStrutRemob(NE) = RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*RCCE(NE)
    ENDDO
 
    IF(Root2ndStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
      Frac2Senes2=AZMAX1(AMIN1(1.0_r8,RCO2MaintDef2ndStruct_Oltd/Root2ndStrutRemob(ielmc)))
    ELSE
      Frac2Senes2=1.0_r8
    ENDIF
  ELSE
    Root2ndStrutRemob(1:NumPlantChemElms) = 0._r8
    Frac2Senes2                           = 0._r8
  ENDIF
  !
  !     SECONDARY ROOT LitrFall CAUSED BY REMOBILIZATION
  !
  !     CSNC,ZSNC,PSNC=literfall C,N,P
  !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
  !
  litrflx2=0._r8
  if(Frac2Senes2.GT.0._r8)then 
    DO NE=1,NumPlantChemElms
      Remobl2ndKill(NE) = RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*Frac2Senes2
      Remobl2ndcycl(NE) = Root2ndStrutRemob(NE)*Frac2Senes2
    ENDDO  
    
    !Of the killed biomass, a fraction is added to litter, while the other is added back to reserve
    D6350: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        dsenecE = AZMAX1(Remobl2ndKill(NE)-Remobl2ndcycl(NE))
        dwoodyE = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_woody_comp)
        dfineE  = PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_fine_comp)

        LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+dwoodyE
        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)  = LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dfineE
        litrflx2(NE)                             = litrflx2(NE)+dwoodyE+dfineE

      ENDDO    
    ENDDO D6350

    ! Remobilization is added to nonstrucal metabolites       
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)+Remobl2ndcycl(NE)
      Root2ndNetGrowthElms(NE)          = Root2ndNetGrowthElms(NE)-Remobl2ndKill(NE)
    ENDDO

    !remobilization shrinks roots
    Root2ndPopExtenz=Root2ndPopExtenz-Frac2Senes2*Root2ndLen_rpvr(N,L,NR,NZ)

  endif
  end associate
  end subroutine DiagSenes2ndRootAxes
!----------------------------------------------------------------------------------------------------
  subroutine GrowRootMycoAxes(I,J,N,L,Lnext,NZ,FoundRootAxesTip,Root1stTipUpdateFlag,TFN6_vr,&
    RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,flag2ndGrowth,CNRTW,CPRTW,fRootGrowPSISense,&
    TotPopuRoot2ndLlenAxes,TotRoot2ndPopuC,litrflx,RCO2flx)
  !
  !Description
  !Grow root axes in laye L  
  implicit none
  integer , intent(in) :: I,J
  INTEGER , INTENT(IN) :: N         !root type id
  integer , intent(in) :: L         !current root layer
  integer , intent(in) :: Lnext     !next root growable layer, soil porosity >0 
  integer , intent(in) :: NZ     !pft id
  real(r8), intent(in) :: TFN6_vr(JZ1)
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  logical,  intent(in) :: flag2ndGrowth(JZ1,MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: CNRTW,CPRTW
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(out) :: TotPopuRoot2ndLlenAxes
  real(r8), intent(out) :: TotRoot2ndPopuC  !secondary root carbon
  logical , intent(inout) :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)  
  logical , intent(inout) :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)  
  real(r8), intent(inout) :: litrflx(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flx
  character(len=*), parameter :: subname='GrowRootMycoAxes'

  real(r8) :: Root2ndStrutRemob(NumPlantChemElms)    
  real(r8) :: DMRespEff                          !respiraiton efficiency for new root growth
  real(r8) :: RespElongWatSens  
  real(r8) :: FRCO2
  real(r8) :: FSNCM
  real(r8) :: FSNCP
  real(r8) :: CNG,CPG
  real(r8) :: Nutstress4GrossResp                           !nutrient stress on root mobilization of nonstrucal C.
  integer  :: NR
  logical  :: FoundPrimaryRootsLayer
  real(r8) :: litrflxt(NumPlantChemElms),litrflx2(NumPlantChemElms)
  real(r8) :: RCO2flxt,RCO2flx2
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: massr1st(NumPlantChemElms),massr1st1(NumPlantChemElms)
  real(r8) :: massr2nd(NumPlantChemElms),massr2nd1(NumPlantChemElms)
  real(r8) :: massnonst(NumPlantChemElms),massnonst1(NumPlantChemElms)
  real(r8) :: massnodul(NumPlantChemElms),massnodul1(NumPlantChemElms)
  real(r8) :: dLext1st,dLext2nd          !turgor pressure allowed 1st and 2nd root elongation [m]
  real(r8) :: RootNsink2nd
  real(r8) :: RTN2X,RTN2Y

  !begin_execution
  associate(                                                                   &
    FracRootElmAllocm            => plt_allom%FracRootElmAllocm               ,& !input  :C woody fraction in root,[-]    
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                   ,& !input  :threshold zero for plang growth calculation, [-]    
    ZERO                         => plt_site%ZERO                             ,& !input  :threshold zero for numerical stability, [-]  
    Root1stLenPP_rpvr            => plt_morph%Root1stLenPP_rpvr               ,& !input  :root layer length primary axes, [m d-2]
    Root2ndRadius_rpvr           => plt_morph%Root2ndRadius_rpvr              ,& !input  :root layer radius secondary axes, [m]
    k_fine_comp                  => pltpar%k_fine_comp                        ,& !input  :fine litter complex id
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft           ,& !input  :plant growth type (vascular, non-vascular),[-]
    CumSoilThickness_vr          => plt_site%CumSoilThickness_vr              ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]    
    NRoot1stTipLay_raxes         => plt_morph%NRoot1stTipLay_raxes            ,& !input  :maximum soil layer number for root axes, [-]
    RootBiomGrosYld_pft          => plt_allom%RootBiomGrosYld_pft             ,& !input  :root growth yield, [g g-1]
    NumAxesPerPrimRoot_pft       => plt_morph%NumAxesPerPrimRoot_pft          ,& !input :primary root axes number, [d-2]            
    PlantPopulation_pft          => plt_site%PlantPopulation_pft              ,& !input  :plant population, [d-2]      
    DLYR3                        => plt_site%DLYR3                            ,& !input  :vertical thickness of soil layer, [m]             
    RootNonstructElmConc_rpvr    => plt_biom%RootNonstructElmConc_rpvr        ,& !input  :root layer nonstructural C concentration, [g g-1]    
    Root2ndLen_rpvr              => plt_morph%Root2ndLen_rpvr                 ,& !input  :root layer length secondary axes, [m d-2]
    RootMyco2ndStrutElms_rpvr    => plt_biom%RootMyco2ndStrutElms_rpvr        ,& !input  :root layer element secondary axes, [g d-2]    
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr        ,& !input  :root layer element primary axes, [g d-2]    
    RootBranchFreq_pft           => plt_morph%RootBranchFreq_pft              ,& !input  :root brancing frequency, [m-1]
    PSIRootTurg_vr               => plt_ew%PSIRootTurg_vr                     ,& !input  :root turgor water potential, [Mpa]
    NumPrimeRootAxes_pft         => plt_morph%NumPrimeRootAxes_pft            ,& !input  :root primary axis number,[-]
    SoilBulkModulus4RootPent_vr  => plt_soilchem%SoilBulkModulus4RootPent_vr  ,& !input  :soil hydraulic resistance, [MPa h m-2]
    NMaxRootBotLayer_pft         => plt_morph%NMaxRootBotLayer_pft            ,& !inoput :maximum soil layer number for all root axes, [-]
    Root1stAxesTipDepz2Surf_pft  => plt_morph%Root1stAxesTipDepz2Surf_pft     ,& !output :plant primary depth relative to column surface, [m]     
    ROOTNLim_rpvr                => plt_biom%ROOTNLim_rpvr                    ,& !output :root N-limitation, 0->1 weaker limitation, [-]     
    ROOTPLim_rpvr                => plt_biom%ROOTPLim_rpvr                     & !output :root P-limitation, 0->1 weaker limitation, [-]                 
  )
  call PrintInfo('beg '//subname)
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     SoilBulkModulus4RootPent_vr,SoilRest4Root2ndPentrate=soil resistance to secondary root penetration (MPa)
!     Root2ndRadius_rpvr=secondary root radius
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     fRootGrowPSISense,RespElongWatSens=growth,respiration function of root water potential
!     PSIRoot_pvr,PSIRootTurg_vr=root total,turgor water potential
!     DMRT=root growth yield
!
!   call SumRootBiome(NZ,mass_inital,massr1st1,massr2nd1,massnonst1,massnodul1)
  !set fine root elongation zone to 2 mm, an intermediate value for fine roots and mycorrhizae
  CALL RootElongationWaterFunc(N,L,NZ,2.e-3_r8,dLext2nd,RespElongWatSens)

  RTN2X    = RootBranchFreq_pft(NZ)*NumAxesPerPrimRoot_pft(NZ)    !fine roots
  RTN2Y    = RootBranchFreq_pft(NZ)*RTN2X                       !root hairs
  dLext2nd = dLext2nd*(RTN2X+RTN2Y)*DLYR3(L)                   !assume fine roots and root hairs are similar in experiencing soil stress
  
  if(N.EQ.ipltroot)then
    !10 mm for primary roots
    CALL RootElongationWaterFunc(N,L,NZ,10.e-3_r8,dLext1st)
    dLext1st=dLext1st*NumAxesPerPrimRoot_pft(NZ)/PlantPopulation_pft(NZ)
  endif  

  !respiration fraction
  DMRespEff         = 1.0_r8-RootBiomGrosYld_pft(NZ)     !1-[gC root/gC nonst] = [gC resp/gC nonst]
  TotPopuRoot2ndLlenAxes = 0._r8;TotRoot2ndPopuC      = 0._r8

  !nutrient limitation factor
  IF(RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CNG  = RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmn,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CNKI)
    CPG  = RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmp,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CPKI)
    Nutstress4GrossResp = AMIN1(CNG,CPG)
  ELSE
    CNG =0._r8
    CPG =0._r8
    Nutstress4GrossResp=1.0_r8
  ENDIF
  ROOTNLim_rpvr(N,L,NZ)=CNG
  ROOTPLim_rpvr(N,L,NZ)=CPG

  !     FOR EACH ROOT AXIS
  !=================================================================================
  !loop through all root axes 
  D5050: DO NR=1,NumPrimeRootAxes_pft(NZ)

    !     SECONDARY ROOT EXTENSION
    !   make sure current layer is no deeper than the lowest root layer
    IF(L.LE.NRoot1stTipLay_raxes(NR,NZ) .AND. .not.FoundRootAxesTip(N,NR))THEN
      RootNsink2nd=1._r8    
      if(N.eq.ipltroot)then
        IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
          RootNsink2nd=Root1stSink_pvr(L,NR)/RootSinkC_vr(N,L)
        endif  
      endif

      !secondary roots can grow in any root layer that is not the deepest
      call Grow2ndRootAxes(I,J,L,NZ,N,NR,Nutstress4GrossResp,CNRTW,CPRTW,dLext2nd,RespElongWatSens,&
        DMRespEff,fRootGrowPSISense,TFN6_vr,RootSinkC_vr,Root2ndSink_pvr,RootNsink2nd,litrflx2,RCO2flx2,Root2ndStrutRemob)

      litrflx = litrflx+litrflx2
      RCO2flx = RCO2flx+RCO2flx2

      TotPopuRoot2ndLlenAxes = TotPopuRoot2ndLlenAxes+Root2ndLen_rpvr(N,L,NR,NZ)
      TotRoot2ndPopuC        = TotRoot2ndPopuC+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
      !
      !make sure it is plant root, and root axis NR is below layer L-1.      
      !
      FoundPrimaryRootsLayer=N.EQ.ipltroot .and. Root1stAxesTipDepz2Surf_pft(NR,NZ).GT.CumSoilThickness_vr(L-1)

      IF(FoundPrimaryRootsLayer)THEN 
        ! plant PRIMARY ROOT EXTENSION
        call Grow1stRootAxes(I,J,N,NR,L,Lnext,NZ,Nutstress4GrossResp,Root1stSink_pvr,&
          Root2ndSink_pvr,RootSinkC_vr,Root2ndStrutRemob,fRootGrowPSISense,TFN6_vr,flag2ndGrowth,dLext1st,RespElongWatSens,&
          DMRespEff,CNRTW,CPRTW,RootNsink2nd,Root1stTipUpdateFlag,FoundRootAxesTip,litrflxt,RCO2flxt)

        litrflx = litrflx+litrflxt
        RCO2flx = RCO2flx+RCO2flxt
      ENDIF

    ENDIF
    NMaxRootBotLayer_pft(NZ)=MAX(NMaxRootBotLayer_pft(NZ),NRoot1stTipLay_raxes(NR,NZ))
    !   call SumRootBiome(NZ,mass_finale,massr1st,massr2nd,massnonst,massnodul)

  ENDDO D5050

  call PrintInfo('end '//subname)
  end associate
  end subroutine GrowRootMycoAxes
!----------------------------------------------------------------------------------------------------
  subroutine GetPlantRoot1stDepz(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  INTEGER :: N,NR,L
  associate(                                                                &
    ZERO                         => plt_site%ZERO                          ,& !input  :threshold zero for numerical stability, [-]  
    NU                           => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    MaxNumRootLays               => plt_site%MaxNumRootLays                ,& !input  :maximum root layer number,[-]    
    NumPrimeRootAxes_pft         => plt_morph%NumPrimeRootAxes_pft         ,& !input  :root primary axis number,[-]    
    Root1stDepz_raxes            => plt_morph%Root1stDepz_raxes            ,& !input  :root layer depth, [m]    
    CumSoilThickness_vr          => plt_site%CumSoilThickness_vr           ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]    
    SoilBulkDensity_vr           => plt_soilchem%SoilBulkDensity_vr        ,& !input  :soil bulk density, [Mg m-3]
    Root1stAxesTipDepz2Surf_pft  => plt_morph%Root1stAxesTipDepz2Surf_pft   & !output :plant primary depth relative to column surface, [m] 
  )  

  N=ipltroot
  DO NR=1,NumPrimeRootAxes_pft(NZ)  
    DO L=NU,MaxNumRootLays
      IF(Root1stDepz_raxes(NR,NZ)>CumSoilThickness_vr(L-1) .AND. Root1stDepz_raxes(NR,NZ)<=CumSoilThickness_vr(L))THEN
        IF(SoilBulkDensity_vr(L).GT.ZERO)THEN
          !in soil
          Root1stAxesTipDepz2Surf_pft(NR,NZ)=Root1stDepz_raxes(NR,NZ)-CumSoilThickness_vr(0)
        ELSE
          !in ponding water
          Root1stAxesTipDepz2Surf_pft(NR,NZ)=Root1stDepz_raxes(NR,NZ)
        ENDIF
        EXIT
      ENDIF
    ENDDO
  ENDDO
  end associate
  end subroutine GetPlantRoot1stDepz
!----------------------------------------------------------------------------------------------------
  subroutine GrowRootElongZone(I,J,N,NR,L,Lnext,NZ,Nutstress4GrossResp,Root1stSink_pvr,&
    Root2ndSink_pvr,RootSinkC_vr,Root2ndStrutRemob,fRootGrowPSISense,TFN6_vr,dLext1st,RespElongWatSens,&
    DMRespEff,CNRTW,CPRTW,RootNsink2nd,litrflxt,RCO2flxt)
  implicit none
  integer,  intent(in) :: I,J,N,NR,NZ
  integer,  intent(in) :: L,Lnext       !current and next layer for root growth
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: dLext1st,RespElongWatSens
  real(r8), intent(in) :: Root2ndStrutRemob(NumPlantChemElms)  
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRespEff                !respiraiton efficiency for new dry root mass creation from nonstructural C
  real(r8), intent(in) :: Nutstress4GrossResp                     !N,P constraint on respiration
  real(r8), intent(in) :: RootNsink2nd
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flxt       !accumulator of CO2 flux
  real(r8) :: FracRoot1stCSinkL
  real(r8) :: GRTWTM  
  real(r8) :: TFRCO2,FRCO2  
  real(r8) :: Rmaint1st_CO2
  real(r8) :: RNonstCO2_OUltd
  real(r8) :: RGrowCO2_Oltd,RGrowCO2_OUltd               !growth respiration unltd,ltd by O2 and unlimited by N,P
  real(r8) :: RNonstCO2_Oltd
  real(r8) :: RCO2XMaint1st_OUltd,RCO2XMaint1st_Oltd
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4GrowC_OUltd
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)  !nonstructural N,P unltd, not ltd by O2 used in growth
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms) !nonstructural N,P unltd, ltd by O2 used in growth
  real(r8) :: RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd
  real(r8) :: FSNCM,FSNCP,CNG,CPG  
  real(r8) :: dsenecE
  real(r8) :: RootMycoNonst4GrowC_Oltd
  real(r8) :: dRootMycoElms(NumPlantChemElms)  !nonstrucal biomass element used for metabolism
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms),TipRadius,RLEN,RootNsink1st
  integer  :: NE,M,LL,LL1
  real(r8) :: dRlenL,RRadius,respscal,ratio
  real(r8) :: Root1stPotExtzPP,Root1stExtenzPP             !pontential/actual  root extension due to nonstrucal C mobilization.  
  real(r8), parameter :: RTDPX1=0.005_r8  !the apical growing zone length
  character(len=*), parameter :: subname='GrowRootElongZone'

  associate(                                                                &
    iPlantPhenolType_pft         => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    fTgrowRootP_vr               => plt_pheno%fTgrowRootP_vr               ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr       => plt_rbgc%RAutoRootO2Limter_rpvr        ,& !input  :O2 constraint to root respiration (0-1), [-]
    MaxNumRootLays               => plt_site%MaxNumRootLays                ,& !input  :maximum root layer number,[-]
    CumSoilThickness_vr          => plt_site%CumSoilThickness_vr           ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    CNRTS_pft                    => plt_allom%CNRTS_pft                    ,& !input  :root N:C ratio x root growth yield, [-]
    CPRTS_pft                    => plt_allom%CPRTS_pft                    ,& !input  :root P:C ratio x root growth yield, [-]
    RootBiomGrosYld_pft          => plt_allom%RootBiomGrosYld_pft          ,& !input  :root growth yield, [g g-1]
    DLYR3                        => plt_site%DLYR3                         ,& !input  :vertical thickness of soil layer, [m]        
    Root1stRadius_pvr            => plt_morph%Root1stRadius_pvr            ,& !input  :root layer radius primary axes, [m]  
    rPCStalk_pft                 => plt_allom%rPCStalk_pft                 ,& !input  :stalk P:C ratio, [g g-1]
    rNCStalk_pft                 => plt_allom%rNCStalk_pft                 ,& !input  :stalk N:C ratio, [gN gC-1]          
    StalkBiomGrowthYld_pft       => plt_allom%StalkBiomGrowthYld_pft       ,& !input  :stalk growth yield, [gC gC-1]    
    Root1stMaxRadius1_pft        => plt_morph%Root1stMaxRadius1_pft        ,& !input  :root radius primary axes, [m]    
    RootMatureAge_pft            => plt_morph%RootMatureAge_pft            ,& !input : Root maturation age, [h]
    PSIRoot_pvr                  => plt_ew%PSIRoot_pvr                     ,& !input :root total water potential, [Mpa]    
    RootAge_rpvr                 => plt_morph%RootAge_rpvr                 ,& !inoput :root age,[h]
    PlantPopulation_pft          => plt_site%PlantPopulation_pft           ,& !input  :plant population, [d-2]
    Root1stSpecLen_pft           => plt_morph%Root1stSpecLen_pft           ,& !input  :specific root length primary axes, [m root gC-1]
    k_fine_comp                  => pltpar%k_fine_comp                     ,& !input  :fine litter complex id    
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr              ,& !input  :root layer length primary axes, [m d-2]
    SeedDepth_pft                => plt_morph%SeedDepth_pft                ,& !input  :seeding depth, [m]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft           ,& !input  :soil layer at planting depth, [-]
    Root1stDepz_raxes            => plt_morph%Root1stDepz_raxes            ,& !input  :root layer depth, [m]
    Root1stMaxRadius_pft         => plt_morph%Root1stMaxRadius_pft         ,& !input  :maximum radius of primary roots, [m]    
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft            ,& !input  :number of main branch,[-]
    FracRootElmAllocm            => plt_allom%FracRootElmAllocm            ,& !input  :C woody fraction in root,[-]    
    Root1stAxesTipDepz2Surf_pft  => plt_morph%Root1stAxesTipDepz2Surf_pft  ,& !output :plant primary depth relative to column surface, [m]     
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft        ,& !input  :plant growth type (vascular, non-vascular),[-]
    GrainFillDowreg_brch         => plt_photo%GrainFillDowreg_brch         ,& !input  :grain fill down-regulation of annual plants, [-]
    RootMaintDef_CO2_pvr         => plt_bgcr%RootMaintDef_CO2_pvr          ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]    
    RootCO2Autor_pvr             => plt_rbgc%RootCO2Autor_pvr              ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr           => plt_rbgc%RootCO2EmisPot_pvr            ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr           => plt_rbgc%RootRespPotent_pvr            ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr        ,& !inoput :root layer nonstructural element, [g d-2]
    RootMyco1stElm_raxs          => plt_biom%RootMyco1stElm_raxs           ,& !inoput :root C primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr    => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !inoput :root layer element secondary axes, [g d-2]
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr     ,& !inoput :root layer element primary axes, [g d-2]
    NRoot1stTipLay_raxes         => plt_morph%NRoot1stTipLay_raxes          & !inoput :soil layer number for deepest root axes, [-]
  )
  call PrintInfo('beg '//subname)
  !
  !     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
  !     ATTRIBUTED TO CURRENT AXIS
  !
  !     Root1stSink_pvr=primary root sink strength
  !     RootSinkC_vr=total root sink strength
  !     FracRoot1stCSinkL=fraction of primary root sink strength in axis
  !

  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot1stCSinkL=Root1stSink_pvr(L,NR)/RootSinkC_vr(N,L)
  ELSE
    FracRoot1stCSinkL=1.0_r8
  ENDIF
  RootNsink1st=1._r8-RootNsink2nd
  !
  !     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
  !     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
  !
  TipRadius = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
  dRlenL    = Root1stAxesTipDepz2Surf_pft(NR,NZ)-CumSoilThickness_vr(L-1)

  !
  !     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
  !     NON-STRUCTURAL C:N:P
  !
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     
  !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !
  !     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
  !     ROOT STRUCTURAL N
  !
  !     Rmaint1st_CO2=root maintenance respiration
  !     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
  !     WTRT1N=primary root N mass
  !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  !     fRootGrowPSISense=growth function of root water potential
  !
  Rmaint1st_CO2 = AZMAX1(RmSpecPlant*RootMyco1stStrutElms_rpvr(ielmn,L,NR,NZ))*TFN6_vr(L)

  IF(is_root_shallow(iPlantRootProfile_pft(NZ)) .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint1st_CO2=Rmaint1st_CO2*fRootGrowPSISense
  ENDIF
  !
  !     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
  !     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
  !
  !     RNonstCO2_OUltd=respiration from non-structural C unlimited by O2
  !     VMXC=rate constant for nonstructural C oxidation in respiration C     
  !     FracRoot1stCSinkL=fraction of primary root sink strength in axis
  !     CPOOL=non-structural C mass
  !     fTgrowRootP_vr=temperature function for root growth
  !     GrainFillDowreg_brch=termination feedback inhibition on C3 CO2
  !     fRootGrowPSISense=growth function of root water potential
  !
  respscal=VMXC*FracRoot1stCSinkL*fTgrowRootP_vr(L,NZ)*Nutstress4GrossResp*GrainFillDowreg_brch(MainBranchNum_pft(NZ),NZ)*fRootGrowPSISense

  RNonstCO2_OUltd=AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))*respscal 

  IF(Root1stAxesTipDepz2Surf_pft(NR,NZ).GE.CumSoilThickness_vr(MaxNumRootLays))THEN
    RNonstCO2_OUltd=AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)
  ENDIF
  !
  !     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
  !
  !     RNonstCO2_Oltd=respiration from non-structural C limited by O2
  !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !     RCO2XMaint1st_OUltd,RCO2XMaint1st_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
  !     RespElongWatSens=respiration function of root water potential
  !
  RNonstCO2_Oltd               = RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RCO2XMaint1st_OUltd          = RNonstCO2_OUltd-Rmaint1st_CO2
  RCO2XMaint1st_Oltd           = RNonstCO2_Oltd-Rmaint1st_CO2
  RGrowCO2_OUltd               = AZMAX1(RCO2XMaint1st_OUltd)*RespElongWatSens
  RGrowCO2_Oltd                = AZMAX1(RCO2XMaint1st_Oltd)*RespElongWatSens
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RCO2XMaint1st_Oltd,0._r8)
  !
  !     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
  !     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
  !
  !     FracRoot1stCSinkL=fraction of secondary root sink strength in axis
  !     ZPOOLR,PPOOLR=non-structural N,P mass in root 
  !     CNRTS_pft,CPRTS_pft=N,P root growth yield = = [gN,P root/gC root]*[gC root /gC nonst]=[gN,P root/gC nonstC]
  !     FNP=growth respiration limited by non-structural N,P
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !
  DMRTR  = DMRespEff     ![gC CO2/gC nonstC]
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*FracRoot1stCSinkL
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*FracRoot1stCSinkL
  FNP    = AMIN1(ZPOOLB/CNRTS_pft(NZ),PPOOLB/CPRTS_pft(NZ))*DMRTR

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
  !     DMRespEff=root C respiration vs nonstructural C consumption
  !     RootMycoNonst4Grow_OUltd(ielmc),RootMycoNonst4Grow_Oltd(ielmc)=root C growth unltd,ltd by O2
  !     DMRT=root growth yield
  !     
  !     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !
  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRespEff
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRespEff
  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)
  !nonstrucal biomass elements for metabolism
  RootMycoNonst4Grow_Oltd(ielmc)  = RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn)  = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
  RootMycoNonst4Grow_Oltd(ielmp)  = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))
  RCO2Nonst4Nassim_OUltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_OUltd(ielmn))
  RCO2Nonst4Nassim_Oltd           = AZMAX1(1.70_r8*RootMycoNonst4Grow_Oltd(ielmn))

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
  dRootMycoElms  = 0._r8
  RCO2T1st_OUltd = AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd+RCO2Nonst4Nassim_OUltd
  RCO2T1st_Oltd  = AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+RCO2Nonst4Nassim_Oltd

  Root1stPotExtzPP = RootMycoNonst4Grow_Oltd(ielmc)/PlantPopulation_pft(NZ)*FracRootElmAllocm(ielmc,k_fine_comp)*Root1stSpecLen_pft(N,NZ)
  Root1stExtenzPP  = Root1stPotExtzPP

  if(lsoilCompaction)then
    if(Root1stPotExtzPP>0._r8)then
      Root1stExtenzPP  = AMIN1(Root1stPotExtzPP,dLext1st)
      if(Root1stPotExtzPP>Root1stExtenzPP)then
        DO NE=1,NumPlantChemElms
          ratio=Root1stExtenzPP/Root1stPotExtzPP
          RootMycoNonst4Grow_Oltd(NE)=RootMycoNonst4Grow_Oltd(NE)*ratio
        ENDDO
      endif  
    else
      Root1stExtenzPP = 0._r8
    endif
  endif

  !dealing with negative growth

  Root1stNetGrowthElms(:)=RootMycoNonst4Grow_Oltd(:)

  call RemobilizePrimeRoots(N,L,NZ,NR,RCO2XMaint1st_Oltd,RCO2XMaint1st_OUltd,dRootMycoElms,&
    RCO2T1st_OUltd,RCO2T1st_Oltd,Root1stNetGrowthElms,litrflxt,RootTipZone=.True.)

  DO NE=1,NumPlantChemElms
    dRootMycoElms(NE)                 = dRootMycoElms(NE)-RootMycoNonst4Grow_Oltd(NE)
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)+dRootMycoElms(NE)
  ENDDO         

  !
  !     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
  !     THROUGH WHICH PRIMARY ROOTS GROW
  !
  !     Root1stLenPP_rpvr=primary root length, assuming primary roots are of the same radius
  !     SeedDepth_pft=seeding depth
  !     FRCO2=fraction of primary root respiration attributed to layer
  !     RootCO2Autor_pvr=total root respiration
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  ! the CO2 from the primary root tip allocated to all layers above the tip, as a temporary fix
  ! for the lack of growth in the maturation zone
  IF(Root1stDepz_raxes(NR,NZ).GT.CumSoilThickness_vr(NGTopRootLayer_pft(NZ)))THEN

    TFRCO2=0._r8;RLEN=0._r8
    LL1=NGTopRootLayer_pft(NZ)
    !if it is a tree, remove mature layers
    if(is_plant_treelike(iPlantRootProfile_pft(NZ)))then 
      !Identify the young layers
      DO LL1=NGTopRootLayer_pft(NZ),NRoot1stTipLay_raxes(NR,NZ)
        IF(RootAge_rpvr(LL1,NR,NZ)<RootMatureAge_pft(NZ))exit
        RLEN=RLEN+Root1stLenPP_rpvr(LL1,NR,NZ)
      ENDDO
    endif
    D5100: DO LL=LL1,NRoot1stTipLay_raxes(NR,NZ)
      !only distribute over young roots
      IF(LL.LT.NRoot1stTipLay_raxes(NR,NZ))THEN
        FRCO2=AMIN1(1.0_r8,Root1stLenPP_rpvr(LL,NR,NZ)/(Root1stDepz_raxes(NR,NZ)-SeedDepth_pft(NZ)-RLEN))
      ELSE
        FRCO2=1.0_r8-TFRCO2
      ENDIF
      TFRCO2                      = TFRCO2+FRCO2
      RootRespPotent_pvr(N,LL,NZ) = RootRespPotent_pvr(N,LL,NZ)+RCO2T1st_OUltd*FRCO2
      RootCO2EmisPot_pvr(N,LL,NZ) = RootCO2EmisPot_pvr(N,LL,NZ)+RCO2T1st_Oltd*FRCO2
      RootCO2Autor_pvr(N,LL,NZ)   = RootCO2Autor_pvr(N,LL,NZ)-RCO2T1st_Oltd*FRCO2
      RCO2flxt                    = RCO2flxt-RCO2T1st_Oltd*FRCO2
    ENDDO D5100
  ELSE
    RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T1st_OUltd
    RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T1st_Oltd
    RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T1st_Oltd
    RCO2flxt                   = RCO2flxt-RCO2T1st_Oltd
  ENDIF
  !
  !     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
  !     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
  !     HAVE DISAPPEARED
  !
  !     RootMycoNonst4Grow_Oltd(ielmc)=primary root C growth ltd by O2
  !     Frac2Senes1=fraction of primary root C to be remobilized
  !     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !     RootMycoNonst4Grow_Oltd(ielmn),RootMycoNonst4Grow_Oltd(ielmp)=nonstructural N,P ltd by O2 used in growth
  !     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
  !     Root2ndLen_rpvr=secondary root length

  !if primary root negative growth is negative, then withdraw secondary roots
  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN          
    call Withdraw2ndRoots(N,NZ,L,NR,Root1stNetGrowthElms,litrflxt)
  ENDIF

  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8 .AND. RootMyco1stElm_raxs(ielmc,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    !primary roots withdraw, note that primary root depth was initialized at seedDepth
    Root1stExtenzPP=Root1stExtenzPP+Root1stNetGrowthElms(ielmc) &
      *(Root1stDepz_raxes(NR,NZ)-SeedDepth_pft(NZ))/RootMyco1stElm_raxs(ielmc,NR,NZ)
  ENDIF

  !the extension should not exceed soil layer thickness
  IF(L.LT.MaxNumRootLays)THEN
    Root1stExtenzPP=AMIN1(DLYR3(Lnext),Root1stExtenzPP)
  ENDIF

  call PrimeRootsExtension(I,J,L,Lnext,N,NR,NZ,FracRoot1stCSinkL,Root1stNetGrowthElms,Root1stExtenzPP)

  call PrintInfo('end '//subname)
  end associate
  end subroutine GrowRootElongZone

!----------------------------------------------------------------------------------------------------
  subroutine Grow1stRootAxes(I,J,N,NR,L,Lnext,NZ,Nutstress4GrossResp,Root1stSink_pvr,&
    Root2ndSink_pvr,RootSinkC_vr,Root2ndStrutRemob,fRootGrowPSISense,TFN6_vr,flag2ndGrowth,dLext1st,RespElongWatSens,&
    DMRespEff,CNRTW,CPRTW,RootNsink2nd,Root1stTipUpdateFlag,FoundRootAxesTip,litrflxt,RCO2flxt)
  !
  !Description:
  !grow primary root axes
  implicit none
  integer,  intent(in) :: I,J,N,NR,NZ
  integer,  intent(in) :: L,Lnext        !current and next layer for root growth
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: dLext1st,RespElongWatSens
  real(r8), intent(in) :: Root2ndStrutRemob(NumPlantChemElms)  
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRespEff              !respiraiton efficiency for fresh root growth
  real(r8), intent(in) :: Nutstress4GrossResp
  real(r8), intent(in) :: RootNsink2nd
  logical, intent(in) :: flag2ndGrowth(JZ1,MaxNumRootAxes)  
  logical,  intent(inout) :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)
  logical,  intent(inout) :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)  
  real(r8), intent(out) :: litrflxt(NumPlantChemElms)
  real(r8), intent(out) :: RCO2flxt
  real(r8) :: Frac2Senes1
  real(r8) :: GRTWTM  
  real(r8) :: TFRCO2,FRCO2  
  real(r8) :: RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd
  real(r8) :: FSNCM,FSNCP,CNG,CPG  
  real(r8) :: dsenecE,SoilResit4RootPentration
  real(r8) :: dRootMycoElms(NumPlantChemElms)  
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: DMTapRTD,dmass(NumPlantChemElms)
  integer  :: NE,M,LL
  logical  :: FoundRootTipLayer
  character(len=*), parameter :: subname='Grow1stRootAxes'

  associate(                                                                &
    iPlantPhenolType_pft         => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    MaxNumRootLays               => plt_site%MaxNumRootLays                ,& !input  :maximum root layer number,[-]
    CumSoilThickness_vr          => plt_site%CumSoilThickness_vr           ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    NumAxesPerPrimRoot_pft       => plt_morph%NumAxesPerPrimRoot_pft       ,& !input  :number of primary root axes per NR, [d-2]        
    Root1stLenPP_rpvr            => plt_morph%Root1stLenPP_rpvr            ,& !input  :root layer length primary axes, [m d-2]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft           ,& !input  :soil layer at planting depth, [-]
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft        ,& !input  :plant growth type (vascular, non-vascular),[-]
    RootMatureAge_pft            => plt_morph%RootMatureAge_pft            ,& !input : Root maturation age, [h]
    RootAge_rpvr                 => plt_morph%RootAge_rpvr                 ,& !inoput :root age,[h]
    Root1stAxesTipDepz2Surf_pft  => plt_morph%Root1stAxesTipDepz2Surf_pft  ,& !output :plant primary depth relative to column surface, [m]     
    RootCO2Autor_pvr             => plt_rbgc%RootCO2Autor_pvr              ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr           => plt_rbgc%RootCO2EmisPot_pvr            ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr           => plt_rbgc%RootRespPotent_pvr            ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr        ,& !inoput :root layer nonstructural element, [g d-2]
    RootMyco1stElm_raxs          => plt_biom%RootMyco1stElm_raxs           ,& !inoput :root chemical element for primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr    => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !inoput :soil layer cheimcal element of secondary roots, [g d-2]
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr     ,& !inoput :soil layer chemical element of primary roots, [g d-2]
    Root1stXNumL_pvr             => plt_morph%Root1stXNumL_pvr             ,& !inoput :soil layer number of primary root axes, [d-2]
    NRoot1stTipLay_raxes         => plt_morph%NRoot1stTipLay_raxes          & !inoput :soil layer number for deepest root axes, [-]
  )
  !   call SumRootBiome(NZ,mass_inital)

  call PrintInfo('beg '//subname)

  litrflxt=0._r8;RCO2flxt=0._r8
  !identify root tip and do root elongation
  !root pass the top surface of layer L, and the primary root referenced by (N,NR) has not been updated  

  FoundRootTipLayer=Root1stAxesTipDepz2Surf_pft(NR,NZ).LE.CumSoilThickness_vr(L) .OR. L.EQ.MaxNumRootLays

  if(.not.Root1stTipUpdateFlag(NR))then
    !diagnose primary root number in layer L
    Root1stXNumL_pvr(L,NZ) = Root1stXNumL_pvr(L,NZ)+NumAxesPerPrimRoot_pft(NZ)
    RootAge_rpvr(L,NR,NZ)  = RootAge_rpvr(L,NR,NZ)+1._r8

    DO NE=1,NumPlantChemElms
      dmass(NE)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMycoNonstElms_rpvr(NE,N,Lnext,NZ)+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ) &
        +RootMyco1stStrutElms_rpvr(NE,Lnext,NR,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)
    ENDDO

    if(FoundRootTipLayer .and. .not. flag2ndGrowth(L,NR))THEN
      !check flat indicating that primary root elongation zone will be updated.
      Root1stTipUpdateFlag(NR)=.true.
      call GrowRootElongZone(I,J,N,NR,L,Lnext,NZ,Nutstress4GrossResp,Root1stSink_pvr,&
        Root2ndSink_pvr,RootSinkC_vr,Root2ndStrutRemob,fRootGrowPSISense,TFN6_vr,dLext1st,RespElongWatSens,&
        DMRespEff,CNRTW,CPRTW,RootNsink2nd,litrflxt,RCO2flxt)
    elseif(flag2ndGrowth(L,NR) .and. Root1stLenPP_rpvr(L,NR,NZ)>0._r8)then
      !behind root tip layer    
      call SecondaryGrowthZone(I,J,N,NR,L,NZ,Nutstress4GrossResp,RootSinkC_vr,Root1stSink_pvr,fRootGrowPSISense,&
        TFN6_vr,RootNsink2nd,litrflxt,RCO2flxt)
    endif

    !consider root withraw at root tip
    IF(L.EQ.NRoot1stTipLay_raxes(NR,NZ))THEN
      call PrimeRootsWithdraw(L,NR,NZ,N,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr)
    ENDIF

    !  REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
    IF(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)       = AZERO(RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))
        RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ) = 0._r8
      ENDDO
    ENDIF

    !  what if nonstructural C is insufficient to meet the negative root mass
    IF(RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ).LT.0.0_r8)THEN
      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)     = RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)
        RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)
        RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = 0._r8
      ENDDO
    ENDIF

    NRoot1stTipLay_raxes(NR,NZ) = MIN(NRoot1stTipLay_raxes(NR,NZ),MaxNumRootLays)
    FoundRootAxesTip(N,NR)      = (L.EQ.NRoot1stTipLay_raxes(NR,NZ))

    DO NE=1,NumPlantChemElms
      dmass(NE)=dmass(NE)-(RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMycoNonstElms_rpvr(NE,N,Lnext,NZ)+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ) &
        +RootMyco1stStrutElms_rpvr(NE,Lnext,NR,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))-litrflxt(NE)
    ENDDO
    dmass(ielmc)=dmass(ielmc)+RCO2flxt
  ENDIF
!  write(666,*)I*1000+J/24.,NR,'L=',L,Lnext,dmass,RCO2flxt
  call PrintInfo('end '//subname)
  end associate
  end subroutine Grow1stRootAxes
!----------------------------------------------------------------------------------------------------

  subroutine SecondaryGrowthZone(I,J,N,NR,L,NZ,Nutstress4GrossResp,RootSinkC_vr,Root1stSink_pvr,fRootGrowPSISense,TFN6_vr,RootNsink2nd,litrflxt,RCO2flxt)
  !
  !Description:
  !Do secondary growth of roots
  implicit none
  integer, intent(in) :: I,J,N,NR,L,NZ
  real(r8), intent(in) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: Nutstress4GrossResp
  real(r8), intent(in) :: RootNsink2nd
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flxt
  real(r8) :: DMRespEff
  real(r8) :: FracRoot1stCSinkL
  real(r8) :: Rmaint1st_CO2
  real(r8) :: RNonstCO2_OUltd,RGrowCO2_Oltd,RGrowCO2_OUltd
  real(r8) :: RNonstCO2_Oltd,RCO2XMaint1st_OUltd,RCO2XMaint1st_Oltd
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4GrowC_Oltd
  real(r8) :: RootMycoNonst4GrowC_OUltd
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)    !nonstrucal biomass converted into structrual biomass with O2 limitation
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms)   !nonstrucal biomass converted into structrual biomass without O2 limitation
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: dRootMycoElms(NumPlantChemElms)              !nonstrucal biomass used for metabolic activity
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd,RootNsink1st,respscal
  real(r8) :: RootThickenWatSens,RespThickenWatSens,dR1stExp,dRPotexp,dRExp,ratio
  integer  :: NE

  character(len=*), parameter :: subname='SecondaryGrowthZone'

  associate(                                                                &
    iPlantPhenolType_pft         => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft            ,& !input  :number of main branch,[-]
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr     ,& !input :root layer element primary axes, [g d-2]
    StalkBiomGrowthYld_pft       => plt_allom%StalkBiomGrowthYld_pft       ,& !input  :stalk growth yield, [gC gC-1]    
    Root1stRadius_pvr            => plt_morph%Root1stRadius_pvr            ,& !input  :root layer radius primary axes, [m]  
    rNCRoot_pft                  => plt_allom%rNCRoot_pft                  ,& !input  :root N:C ratio, [gN gC-1]
    Root1stRadius_rpvr           => plt_morph%Root1stRadius_rpvr           ,& !input: root layer radius for each primary axes,  [m]    
    rPCStalk_pft                 => plt_allom%rPCStalk_pft                 ,& !input  :stalk P:C ratio, [g g-1]
    CoarseRootVolPerMassC_pft    => plt_morph%CoarseRootVolPerMassC_pft    ,& !input  :coarse root volume:mass ratio, [m3 gC-1]    
    rNCStalk_pft                 => plt_allom%rNCStalk_pft                 ,& !input  :stalk N:C ratio, [gN gC-1]          
    GrainFillDowreg_brch         => plt_photo%GrainFillDowreg_brch         ,& !input  :grain fill down-regulation of annual plants, [-]
    PlantPopulation_pft          => plt_site%PlantPopulation_pft           ,& !input  :plant population, [d-2]    
    Root1stLenPP_rpvr            => plt_morph%Root1stLenPP_rpvr            ,& !input :root layer length primary axes per plant, [m d-2]    
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    fTgrowRootP_vr               => plt_pheno%fTgrowRootP_vr               ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr       => plt_rbgc%RAutoRootO2Limter_rpvr        ,& !input  :O2 constraint to root respiration (0-1), [-]
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr        ,& !inoput :root layer nonstructural element, [g d-2]
    RootCO2EmisPot_pvr           => plt_rbgc%RootCO2EmisPot_pvr            ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootMaintDef_CO2_pvr         => plt_bgcr%RootMaintDef_CO2_pvr          ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]    
    RootCO2Autor_pvr             => plt_rbgc%RootCO2Autor_pvr              ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootMyco1stElm_raxs          => plt_biom%RootMyco1stElm_raxs           ,& !inoput :root C primary axes, [g d-2]
    RootCRRadius0_rpvr           => plt_morph%RootCRRadius0_rpvr           ,& !inoput :initial radius of roots that may undergo secondary growth, [m]
    RootRespPotent_pvr           => plt_rbgc%RootRespPotent_pvr             & !inoput :root respiration unconstrained by O2, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)
    !if the plant has deep roots, aka taproot 
  !grow the maturation zone
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot1stCSinkL=Root1stSink_pvr(L,NR)/RootSinkC_vr(N,L)
  ELSE
    FracRoot1stCSinkL=1.0_r8
  ENDIF
  RootNsink1st=1._r8-RootNsink2nd

  !Do maintenance respiration 
  Rmaint1st_CO2 = AZMAX1(RmSpecPlant*RootMyco1stStrutElms_rpvr(ielmn,L,NR,NZ))*TFN6_vr(L)*0.2_r8

  if(iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint1st_CO2=Rmaint1st_CO2*fRootGrowPSISense
  endif

  call RootThickenWaterFunc(N,L,NZ,Root1stRadius_rpvr(L,NR,NZ),RespThickenWatSens,dR1stExp)

  !mobilize storage carbon and 
  respscal=VMXC*FracRoot1stCSinkL*fTgrowRootP_vr(L,NZ)*Nutstress4GrossResp*GrainFillDowreg_brch(MainBranchNum_pft(NZ),NZ) &
    *fRootGrowPSISense*rNCStalk_pft(NZ)/rNCRoot_pft(NZ)
  RNonstCO2_OUltd=AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))*respscal*0.2_r8 

  !apply oxygen dependence
  RNonstCO2_Oltd               = RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  !diagnose growth respiraiton
  RCO2XMaint1st_OUltd          = RNonstCO2_OUltd-Rmaint1st_CO2
  RCO2XMaint1st_Oltd           = RNonstCO2_Oltd-Rmaint1st_CO2
  RGrowCO2_OUltd               = AZMAX1(RCO2XMaint1st_OUltd)
  RGrowCO2_Oltd                = AZMAX1(RCO2XMaint1st_Oltd)
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RCO2XMaint1st_Oltd,0._r8)

  !apply nutrient limitation to growth respiraiton
  DMRespEff  = 1._r8-StalkBiomGrowthYld_pft(NZ)       ![gC CO2/gC nonst]
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*FracRoot1stCSinkL
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*FracRoot1stCSinkL
  FNP    = AMIN1(ZPOOLB/rNCStalk_pft(NZ),PPOOLB/rPCStalk_pft(NZ))*DMRespEff/StalkBiomGrowthYld_pft(NZ)
  
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

  !the biomass yield considered to be simlar for stalk 
  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRespEff
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRespEff
  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*StalkBiomGrowthYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc))*rNCStalk_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc))*rPCStalk_pft(NZ)

  RootMycoNonst4Grow_Oltd(ielmc) = RootMycoNonst4GrowC_Oltd*StalkBiomGrowthYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn) = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*rNCStalk_pft(NZ)))
  RootMycoNonst4Grow_Oltd(ielmp) = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*rPCStalk_pft(NZ)))  

  dRPotexp = RootMycoNonst4Grow_Oltd(ielmc)*CoarseRootVolPerMassC_pft(NZ)/(TwoPiCON*Root1stRadius_rpvr(L,NR,NZ)*Root1stLenPP_rpvr(L,NR,NZ))
  dRExp    = dRPotexp
  if(lsoilCompaction)then
    if(dRPotexp>0._r8)then
      dRExp    = AMIN1(dRPotexp,dR1stExp)
      if(dRPotexp>dRExp)then
        ratio=safe_adb(dRExp,dRPotexp)
        DO NE=1,NumPlantChemElms
          RootMycoNonst4Grow_Oltd(NE)=RootMycoNonst4Grow_Oltd(NE)*ratio
        ENDDO
      endif  
    else
      dRExp=0._r8  
    endif
  endif
  Root1stNetGrowthElms(:)  = RootMycoNonst4Grow_Oltd(:)

  dRootMycoElms  = 0._r8
  RCO2T1st_OUltd = AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd
  RCO2T1st_Oltd  = AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd

  call RemobilizePrimeRoots(N,L,NZ,NR,RCO2XMaint1st_Oltd,RCO2XMaint1st_OUltd,dRootMycoElms,&
    RCO2T1st_OUltd,RCO2T1st_Oltd,Root1stNetGrowthElms,litrflxt,RootTipZone=.false.)

  dRootMycoElms(ielmc)=dRootMycoElms(ielmc)-RCO2T1st_Oltd

  DO NE=1,NumPlantChemElms
    dRootMycoElms(NE)                 = dRootMycoElms(NE)-RootMycoNonst4Grow_Oltd(NE)
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)+dRootMycoElms(NE)
  ENDDO 

  RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T1st_OUltd
  RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T1st_Oltd
  RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T1st_Oltd
  RCO2flxt                   = RCO2flxt-RCO2T1st_Oltd

  !if primary root negative growth is negative, then withdraw secondary roots
  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN          
    call Withdraw2ndRoots(N,NZ,L,NR,Root1stNetGrowthElms,litrflxt)
  ENDIF

  !apply the structrual growth
  DO NE=1,NumPlantChemElms
    RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+Root1stNetGrowthElms(NE)
    RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)+Root1stNetGrowthElms(NE)
  ENDDO

  Root1stRadius_rpvr(L,NR,NZ)=sqrt(RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)*CoarseRootVolPerMassC_pft(NZ) &
    /(PiCON*Root1stLenPP_rpvr(L,NR,NZ)*PlantPopulation_pft(NZ)))
  call PrintInfo('end '//subname)  
  end associate
  end subroutine SecondaryGrowthZone

!----------------------------------------------------------------------------------------------------
  subroutine RootElongationWaterFunc(N,L,NZ,Lelonz,dLext,RespElongWatSens)
  !
  !Description:
  !water dependence for root elongation
  implicit none
  integer, intent(in) :: N,L,NZ
  real(r8), intent(in) :: Lelonz                          !length of the elongation zone. [m]       
  real(r8), intent(out) :: dLext                          !pressure-allowed root penetration/elongation [h m-1]
  real(r8), optional, intent(out) :: RespElongWatSens     !respiration function of root water potential, [-]

  !--- local variables ---
  real(r8), parameter :: phi_lmin = 0.1_r8                !minimum root elastic modulus normalizer [MPa-1 h], i.e. F' in Rickman et al. 1992, Grant, 1992
  real(r8), parameter :: phi_lmax = 1._r8                 !maximum root elastic modulus normalizer [MPa-1 h]
  real(r8), parameter :: Ksens=0.8_r8                     !Haf saturation constant, [MPa]

  real(r8)  :: RootElongWatSens                           !water function for root elongation, [-]
  real(r8) :: V0,V0p                                      !dispaceable soil volume due to root growth [m3]  
  real(r8) :: rv0,Radius_Impact                           !impact radius, [m]
  real(r8), parameter :: r0=0.05e-2_r8                    !reference radius [m]
  real(r8), parameter :: Vr0=PICON*1.e-2_r8*r0**2         !reference root volume, [m3] 
  real(r8) :: phi_l,SoilResit4RootPentration
  character(len=*), parameter :: subname='RootElongationWaterFunc'

  associate(                                                                    &
    PSIRootTurg_vr               => plt_ew%PSIRootTurg_vr                      ,& !input  :root turgor water potential, [Mpa]
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    SoilBulkModulus4RootPent_vr  => plt_soilchem%SoilBulkModulus4RootPent_vr    & !input  :elastic modulus of the undisturbed soil, [MPa]
  )
!     SoilBulkModulus4RootPent_vr,SoilResit4PrimRootPentration=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
  !eq.(3) in Rickman et al (1992), Calculating daily root length density profiles by applying elastic theory to agricultural soils
  
  phi_l=phi_lmin+(phi_lmax-phi_lmin)/(1._r8+(SoilBulkModulus4RootPent_vr(L)/Ksens)**2)
  !eq. (10) in Rickman et al. (1992)
  !V0 was obtained by regression, assuming r0=0.05 cm=0.5mm and elongation zone 1cm.
  V0=(0.0327_r8+0.0055_r8*SoilBulkModulus4RootPent_vr(L))*1.e-4_r8

  !!aka SoilResit4RootPentration=E'*dt/V0 in eq. (4) of Rickman et al. (1992), [h m-1]
  SoilResit4RootPentration = phi_l*SoilBulkModulus4RootPent_vr(L)*Vr0/V0

  !The lockhart equation for elongation
  RootElongWatSens = AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-TurgPSIMin4OrganExtens))  
  if(present(RespElongWatSens))RespElongWatSens = fRespWatSens(RootElongWatSens,iPlantRootProfile_pft(NZ))
  dLext            = RootElongWatSens*phi_l*Lelonz/(1._r8+SoilResit4RootPentration)
  end associate    
  END subroutine RootElongationWaterFunc
!----------------------------------------------------------------------------------------------------
  subroutine RootThickenWaterFunc(N,L,NZ,RootRadius,RespThickenWatSens,dR1stExp)

  implicit none
  integer, intent(in) :: N,L,NZ
  real(r8), intent(in)  :: RootRadius                     !root radius [m]
  real(r8), intent(out) :: RespThickenWatSens             !respiration function of root water potential, [-]
  real(r8), intent(out) :: dR1stExp                       !stress allowed radial expansion, [m]

  real(r8), parameter :: phi_Rmin =0.1_r8      !the baseline extensibility in hydroponics, which keeps root thin [MPa-1 h-1]
  real(r8), parameter :: phi_Rmax =0.5_r8      !The maximum extensibility in hard soil, allows thickening, [MPa-1 h-1]
  real(r8), parameter :: KSensR = 0.1_r8       !The "Sensory Threshold." The pressure at which the plant starts thickening [MPa]
  real(r8) :: phi_R,normSoilStress
  real(r8) :: RootThickenWatSens
  associate(                                                                     &
    PSIRootTurg_vr                => plt_ew%PSIRootTurg_vr                      ,& !input  :root turgor water potential, [Mpa]
    iPlantRootProfile_pft         => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    SoilModulus4RootRadialexp_vr  => plt_soilchem%SoilModulus4RootRadialexp_vr   & !input  :soil modulus for root radial expansion, [MPa]    
  )

  normSoilStress     = (SoilModulus4RootRadialexp_vr(L)/KSensR)**2
  phi_R              = phi_Rmin+(phi_Rmax-phi_Rmin)*(normSoilStress/(1._r8+normSoilStress))
  RootThickenWatSens = AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-TurgPSIMin4OrganExtens-SoilModulus4RootRadialexp_vr(L)))
  RespThickenWatSens = fRespWatSens(RootThickenWatSens,iPlantRootProfile_pft(NZ))
  dR1stExp           = phi_R*RootRadius*RootThickenWatSens
  end associate
  end subroutine RootThickenWaterFunc
!----------------------------------------------------------------------------------------------------
  subroutine Withdraw2ndRoots(N,NZ,L,NR,RootNetGrowthElms,litrflxt)
  !
  !Description
  !Due to negative growth of primary roots, secondary roots are killed to make up the metabolic deficit
  implicit none
  integer, intent(in) :: N        !root type, plant/myco
  integer, intent(in) :: NZ,L,NR
  real(r8), intent(inout) :: RootNetGrowthElms(NumPlantChemElms)
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)

  integer :: LX,LL,NE,M
  real(r8) :: GRTWTM     !negative primary root C growth
  real(r8) :: FSNCM,FSNCP
  real(r8) :: dRootMyco2ndst2Litr(NumPlantChemElms)
  real(r8) :: dRootMycoNonst2Litr(NumPlantChemElms)
  real(r8) :: dsenecE,dwoodyE,dfineE1,dfineE2
  character(len=*), parameter :: subname='Withdraw2ndRoots'


  associate(                                                          &
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]  
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr    ,& !input  :root layer structural C, [gC d-2]
    icwood                    => pltpar%icwood                       ,& !input  :group id of coarse woody litter
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm         ,& !input  :C woody fraction in root,[-]
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    inonstruct                => pltpar%inonstruct                   ,& !input  :group id of plant nonstructural litter
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !inoput :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr            & !inoput :root layer length secondary axes, [m d-2]
  )
  call PrintInfo('beg '//subname)
  LX=MAX(1,L-1)        

  D5105: DO LL=L,LX,-1
    !make a copy of the total negative growth 
    GRTWTM=RootNetGrowthElms(ielmc)
    
    !remove secondary roots and offset some of the negative growth
    D5106: DO NE=1,NumPlantChemElms
      IF(RootNetGrowthElms(NE).LT.0.0_r8)THEN
        IF(RootNetGrowthElms(NE).GT.-RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ))THEN
          if(NE==ielmc)then
            Root2ndLen_rpvr(N,LL,NR,NZ)=Root2ndLen_rpvr(N,LL,NR,NZ)*(1._r8+RootNetGrowthElms(NE)/RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ))
          endif   
          RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ) = RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)+RootNetGrowthElms(NE)
          RootNetGrowthElms(NE)                    = 0._r8
        ELSE
          if(NE==ielmc)Root2ndLen_rpvr(N,LL,NR,NZ)=0._r8
          RootNetGrowthElms(NE)                    = RootNetGrowthElms(NE)+RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)
          RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ) = 0._r8
        ENDIF
      ENDIF
    ENDDO D5106
    !
    !     CONCURRENT LOSS OF MYCORRHIZAE AND NODULES WITH LOSS
    !     OF SECONDARY ROOTS
    !
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
    IF(GRTWTM.LT.0.0_r8 .and. Myco_pft(NZ)>1)THEN
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

      dRootMyco2ndst2Litr = 0._r8
      dRootMycoNonst2Litr = 0._r8

      D6450: DO M=1,jsken
        D64511: DO NE=1,NumPlantChemElms
          dsenecE = FSNCM*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ))
          dwoodyE = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_woody_comp)
          dfineE1 = PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_fine_comp)

          LitrfallElms_pvr(NE,M,k_woody_comp,LL,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,LL,NZ)+dwoodyE
          LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ)  = LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ)+dfineE1
          dRootMyco2ndst2Litr(NE)                   = dRootMyco2ndst2Litr(NE)+dwoodyE+dfineE1

          litrflxt(NE)=litrflxt(NE)+dwoodyE+dfineE1

          dfineE2=PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*FSNCP*(RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ))

          LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ)+dfineE2
          dRootMycoNonst2Litr(NE)                  = dRootMycoNonst2Litr(NE)+dfineE2

          litrflxt(NE)=litrflxt(NE)+dfineE2              
        ENDDO D64511   
      ENDDO D6450

      DO NE=1,NumPlantChemElms  
        RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ) = AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ)-dRootMyco2ndst2Litr(NE))
        RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ)       = AZMAX1(RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ)-dRootMycoNonst2Litr(NE))        
      ENDDO
      Root2ndLen_rpvr(imycorrhz,LL,NR,NZ)=AZMAX1(Root2ndLen_rpvr(imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)
    ENDIF
  ENDDO D5105
  call PrintInfo('end '//subname)
  end associate
  end subroutine Withdraw2ndRoots

!----------------------------------------------------------------------------------------------------
  subroutine RemobilizePrimeRoots(N,L,NZ,NR,RCO2XMaint1st_Oltd,RCO2XMaint1st_OUltd,dRootMycoElms,&
    RCO2T1st_OUltd,RCO2T1st_Oltd,Root1stNetGrowthElms,litrflxt,RootTipZone)

  implicit none
  integer , intent(in) :: N,L,NZ,NR
  real(r8), intent(in) :: RCO2XMaint1st_Oltd                         !diagnostic for maintenance deficit with O2-limitation
  real(r8), intent(in) :: RCO2XMaint1st_OUltd                        !diagnostic for maintenance deficit without O2-limitation
  real(r8), intent(inout) :: RCO2T1st_Oltd                           !CO2 flux from 1st roots, with O-limitation
  real(r8), intent(inout) :: RCO2T1st_OUltd                          !CO2 flux from 1st roots, without O-limitation
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: dRootMycoElms(NumPlantChemElms)
  real(r8), intent(inout) :: Root1stNetGrowthElms(NumPlantChemElms)
  logical,  intent(in)    :: RootTipZone                              !dealing with root tip zone or not?
  integer :: M,NE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: Root1stStrutRemob(NumPlantChemElms)
  real(r8) :: Frac2Senes1,dlitrfall
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: RCO2MaintDef1stStruct_Oltd   !maintenance deficit met by remobilization under O2 limitation
  real(r8) :: RCO2MaintDef1stStruct_OUltd  !maintenance deficit met by remobilization without O2 limitation
  real(r8) :: dsenecE,RCCE(NumPlantChemElms)
  real(r8) :: Root1stcylc(NumPlantChemElms),Root1stKill(NumPlantChemElms)
  CHARACTER(LEN=*), parameter :: subname='RemobilizePrimeRoots'

! begin_execution
  associate(                                                           &
    RootNonstructElmConc_rpvr  => plt_biom%RootNonstructElmConc_rpvr   ,& !input  :root layer nonstructural C concentration, [g g-1]
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    icwood                     => pltpar%icwood                        ,& !input  :group id of coarse woody litter
    istalk                     => pltpar%istalk                        ,& !input  :group id of coarse woody litter
    k_woody_comp               => pltpar%k_woody_comp                  ,& !input  :woody litter complex id
    k_fine_comp                => pltpar%k_fine_comp                   ,& !input  :fine litter complex id
    iroot                      => pltpar%iroot                         ,& !input  :group id of plant root litter
    RootMyco1stStrutElms_rpvr  => plt_biom%RootMyco1stStrutElms_rpvr   ,& !input  :root layer element primary axes, [g d-2]
    FracRootElmAllocm          => plt_allom%FracRootElmAllocm          ,& !input  :root element woody/fine fraction in biomass,[-]
    FracWoodStalkElmAlloc2Litr => plt_allom%FracWoodStalkElmAlloc2Litr ,& !input  :stalk element woody/fine fraction in biomass,[-]
    PlantElmAllocMat4Litr      => plt_soilchem%PlantElmAllocMat4Litr   ,& !input  :litter kinetic fraction, [-]
    RAutoRootO2Limter_rpvr     => plt_rbgc%RAutoRootO2Limter_rpvr      ,& !input  :O2 constraint to root respiration (0-1), [-]
    iPlantCalendar_brch        => plt_pheno%iPlantCalendar_brch        ,& !input  :plant growth stage, [-]
    MainBranchNum_pft          => plt_morph%MainBranchNum_pft          ,& !input  :number of main branch,[-]
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr             & !inoput :plant LitrFall element, [g d-2 h-1]
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
  call PrintInfo('beg '//subname)
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 &
    .AND. RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmn,N,L,NZ), &
      RootNonstructElmConc_rpvr(ielmn,N,L,NZ)+RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CNKI), &
      safe_adb(RootNonstructElmConc_rpvr(ielmp,N,L,NZ),RootNonstructElmConc_rpvr(ielmp,N,L,NZ) &
        +RootNonstructElmConc_rpvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC=0._r8
    CNC=0._r8
    CPC=0._r8
  ENDIF
  RCCC        = RCCZR+CCC*RCCYR
  RCCN        = CNC*RCCXR
  RCCP        = CPC*RCCQR
  RCCE(ielmc) = RCCC
  RCCE(ielmn) = (RCCN+(1.0_r8-RCCN)*RCCC)
  RCCE(ielmp) = (RCCP+(1.0_r8-RCCP)*RCCC)
  !
  !     RECOVERY OF REMOBILIZABLE N,P DURING PRIMARY ROOT REMOBILIZATION
  !     DEPENDS ON ROOT NON-STRUCTURAL C:N:P
  !
  !     RCO2XMaint1st_OUltd,RCO2XMaint1st_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
  !     RCO2MaintDef1stStruct_OUltd,RCO2MaintDef1stStruct_Oltd=excess maintenance respiration unltd,ltd by O2
  !     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !     Root1stStrutRemob(ielmc),Root1stStrutRemob(ielmn),Root1stStrutRemob(ielmp)=remobilization of C,N,P from senescing root
  !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !     Frac2Senes1=fraction of primary root C to be remobilized
  !
  !insufficient energy for maintenance
  IF(-RCO2XMaint1st_OUltd.GT.0.0_r8)THEN
    IF(-RCO2XMaint1st_OUltd.LT.RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)*RCCC)THEN
      RCO2MaintDef1stStruct_OUltd=-RCO2XMaint1st_OUltd
    ELSE
      RCO2MaintDef1stStruct_OUltd=AZMAX1(RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)*RCCC)
    ENDIF
  ELSE
    RCO2MaintDef1stStruct_OUltd=0._r8
  ENDIF

  IF(-RCO2XMaint1st_Oltd.GT.0.0_r8)THEN
    !maintenance deficit is less than remobilizable C.
    IF(-RCO2XMaint1st_Oltd.LT.RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)*RCCC)THEN
      RCO2MaintDef1stStruct_Oltd=-RCO2XMaint1st_Oltd
    ELSE
      RCO2MaintDef1stStruct_Oltd=AZMAX1(RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)*RCCC)*RAutoRootO2Limter_rpvr(N,L,NZ)
    ENDIF
  ELSE
    RCO2MaintDef1stStruct_Oltd=0._r8
  ENDIF

  IF(RCO2MaintDef1stStruct_Oltd.GT.0.0_r8 .AND. RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    !remobilization upon starvation-induced root retreat
    DO NE=1,NumPlantChemElms
      Root1stStrutRemob(NE) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)*RCCE(NE)
    ENDDO

    !maintenance deficit is paid by remobilization
    IF(Root1stStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
      Frac2Senes1=AZMAX1(AMIN1(1.0_r8,RCO2MaintDef1stStruct_Oltd/Root1stStrutRemob(ielmc)))
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
  !
  if(Frac2Senes1>0._r8)then
    DO NE=1,NumPlantChemElms          
      Root1stKill(NE) = Frac2Senes1*RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)  
      Root1stcylc(NE) = Frac2Senes1*Root1stStrutRemob(NE)
    ENDDO  
    if(RootTipZone)then
      D6355: DO M=1,jsken
        DO NE=1,NumPlantChemElms          
          dsenecE      = AZMAX1(Root1stKill(NE)-Root1stcylc(NE))
          dlitrfall    = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_woody_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall

          LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) +dlitrfall
          
          dlitrfall    = PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_fine_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall

          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dlitrfall
        ENDDO
      ENDDO D6355
    else
      D6356: DO M=1,jsken
        DO NE=1,NumPlantChemElms          
          dsenecE      = AZMAX1(Root1stKill(NE)-Root1stcylc(NE))
          dlitrfall    = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*dsenecE*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall

          LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+dlitrfall

          dlitrfall    = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*dsenecE*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall

          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) +dlitrfall
        ENDDO
      ENDDO D6356
    endif
    !remobilization reduces growth

    DO NE=1,NumPlantChemElms
      Root1stNetGrowthElms(NE)=Root1stNetGrowthElms(NE)-Root1stKill(NE)
    ENDDO
  endif
  !
  !     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
  !
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
  !     Rmaint1st_CO2=root maintenance respiration
  !     RNonstCO2_Oltd=respiration from non-structural C limited by O2
  !     RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn ltd by O2
  !     RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !     RCO2MaintDef1stStruct_Oltd=excess maintenance respiration ltd by O2
  !     Frac2Senes1=fraction of primary root C to be remobilized
  !     Root2ndStrutRemob(:)=remobilization of C,N,P from senescing root
  !     RootMycoNonst4Grow_Oltd(ielmn),RootMycoNonst4Grow_Oltd(ielmp)=nonstructural N,P ltd by O2 used in growth
  ! add CO2 respiraiton from remobilization  
  RCO2T1st_Oltd  = RCO2T1st_Oltd+RCO2MaintDef1stStruct_Oltd
  RCO2T1st_OUltd = RCO2T1st_OUltd+RCO2MaintDef1stStruct_OUltd

  DO NE=1,NumPlantChemElms
    dRootMycoElms(NE)=dRootMycoElms(NE)+Root1stcylc(NE) 
  ENDDO 
  dRootMycoElms(ielmc)=dRootMycoElms(ielmc)-RCO2T1st_Oltd

  call PrintInfo('end '//subname)
  end associate
  end subroutine RemobilizePrimeRoots

!----------------------------------------------------------------------------------------------------
  subroutine PrimeRootsExtension(I,J,L,Lnext,N,NR,NZ,FracRoot1stCSinkL,RootNetGrowthElms,Root1stExtenzPP)
  !
  !Description
  !Grow primary roots at the tip

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: L,Lnext  !current and next soil layer to support root growth
  integer, intent(in) :: N    !root/myco indicator
  integer, intent(in) :: NR   !root axis 
  integer, intent(in) :: NZ   !pft number
  real(r8),intent(in) :: FracRoot1stCSinkL
  real(r8), intent(in) :: RootNetGrowthElms(NumPlantChemElms)  !whole population growth for primary roots
  real(r8), intent(in) :: Root1stExtenzPP       !per plant primary root extension
  real(r8) :: FGROL,FGROZ                 !fraction of root extension in current and next lower soil layer
  real(r8) :: FGROLA
  integer :: NE
  REAL(R8) :: XFRE(NumPlantChemElms)
  character(len=*), parameter :: subname='PrimeRootsExtension'

! begin_execution
  associate(                                                          &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    rProteinC2RootN_pft       => plt_allom%rProteinC2RootN_pft       ,& !input  :protein C to root N ratio in root biomass, [-]
    rProteinC2RootP_pft       => plt_allom%rProteinC2RootP_pft       ,& !input  :protein C to root P ratio in root biomass, [-]
    CumSoilThickness_vr       => plt_site%CumSoilThickness_vr        ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    DLYR3                     => plt_site%DLYR3                      ,& !input  :vertical thickness of soil layer, [m]
    MaxNumRootLays            => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft        ,& !input  :soil layer at planting depth, [-]
    SeedDepth_pft             => plt_morph%SeedDepth_pft             ,& !input  :seeding depth, [m]
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr              ,& !input  :root layer volume water, [m2 d-2]    
    Root1stMaxRadius1_pft     => plt_morph%Root1stMaxRadius1_pft     ,& !input  :maximum radius of primary roots, [m]    
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft      ,& !input  :maximum radius of primary roots, [m]    
    RootCRRadius0_rpvr        => plt_morph%RootCRRadius0_rpvr        ,& !inoput: initial radius of roots that may undergo secondary growth, [m]    
    Root1stRadius_rpvr        => plt_morph%Root1stRadius_rpvr        ,& !inoput: root layer radius for each primary axes,  [m]
    RootAge_rpvr              => plt_morph%RootAge_rpvr              ,& !inoput :root age,[h]    
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !inoput :root layer element primary axes, [g d-2]
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs        ,& !inoput :root C primary axes, [g d-2]
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr           ,& !inoput :root layer protein C, [gC d-2]
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr         ,& !inoput :root layer C, [gC d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr               ,& !inoput :root turgor water potential, [Mpa]
    PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr               ,& !inoput :root osmotic water potential, [Mpa]
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr                  ,& !inoput :root total water potential, [Mpa]
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr         ,& !inoput :root layer radius primary axes, [m]
    Root1stLenPP_rpvr         => plt_morph%Root1stLenPP_rpvr         ,& !inoput :root layer length primary axes per plant, [m d-2]
    Root1stDepz_raxes         => plt_morph%Root1stDepz_raxes         ,& !inoput :root layer depth, [m]
    NRoot1stTipLay_raxes      => plt_morph%NRoot1stTipLay_raxes       & !output :soil layer number for deepest root axes, [-]
  )
  call PrintInfo('beg '//subname)
  !  PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
  !
  !  Root1stExtenzPP=primary root length extension
  !  RootMycoNonst4Grow_Oltd=primary root C growth ltd by O2
  !  PP=PFT population
  !  FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
  !  RootNetGrowthElms(:)=net primary root C,N,P growth
  !  RTDP1=primary root depth from soil surface
  !  SeedDepth_pft=seeding depth
  !  Frac2Senes1=fraction of primary root C to be remobilized
  !  Root1stLenPP_rpvr=primary root length
  !  RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !  DLYR=soil layer thickness
  !  only the non-woody part is contributing to the elongation
  ! elongation due to positive growth of non-woody part

  !
  !     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
  !     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
  !     BOUNDARY OF CURRENT LAYER
  !
  !     Root1stExtenzPP=primary root length extension
  !
  ! Question, 03/29/2024, Jinyun Tang: needs double check the following calculation 
  ! if FGROL < 1.0, then the extension is all in current layer, meaning FGROZ=0.0
  !  
  IF(Root1stExtenzPP.GT.ZERO4Groth_pft(NZ) .AND. L.LT.MaxNumRootLays)THEN
    FGROL=AZMAX1(AMIN1(1.0_r8,(CumSoilThickness_vr(L)-Root1stDepz_raxes(NR,NZ))/Root1stExtenzPP))
    IF(FGROL.LT.1.0_r8)FGROL=0._r8
    FGROZ=AZMAX1(1.0_r8-FGROL)
  ELSE
    FGROL = 1.0_r8
    FGROZ = 0._r8
  ENDIF

  FGROLA=AMIN1(1._r8,Root1stExtenzPP*FGROL/DLYR3(L))

  !
  !     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
  !     AND AXIS NUMBER
  !
  !     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !     RootNetGrowthElms(:)=net root C,N,P growth
  !     Root1stExtenzPP=primary root length extension
  !     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
  !     FGROL,FGROZ=fraction of Root1stExtenzPP in current,next lower soil layer
  !     RootProteinC_pvr=total root protein C mass
  !     CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
  !     Root1stLenPP_rpvr=primary root length per plant
  !
  Root1stDepz_raxes(NR,NZ) = Root1stDepz_raxes(NR,NZ)+Root1stExtenzPP
  RootAge_rpvr(L,NR,NZ)    = RootAge_rpvr(L,NR,NZ)*(1._r8-FGROLA)+FGROLA

  DO NE=1,NumPlantChemElms
    RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+RootNetGrowthElms(NE)
    RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)+RootNetGrowthElms(NE)*FGROL
  ENDDO
 
  RootCRRadius0_rpvr(L,NR,NZ) = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
  Root1stRadius_rpvr(L,NR,NZ) = RootCRRadius0_rpvr(L,NR,NZ)
  RootProteinC_pvr(N,L,NZ)    = RootProteinC_pvr(N,L,NZ)+ &
    AMIN1(rProteinC2RootN_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmn,L,NR,NZ),&
       rProteinC2RootP_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmp,L,NR,NZ))
         
  Root1stLenPP_rpvr(L,NR,NZ)=Root1stLenPP_rpvr(L,NR,NZ)+Root1stExtenzPP*FGROL
  !
  !     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
  !     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
  !     OF CURRENT SOIL LAYER
  !
  !     FGROZ=fraction of Root1stExtenzPP into next lower soil layer
  !     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
  !     RootNetGrowthElms(:)=net root C,N,P growth
  !     RootProteinC_pvr=total root protein C mass
  !     CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
  !     WTRTD=root C mass
  !     Root1stLenPP_rpvr=primary root length per plant
  !     Root1stExtenzPP=primary root length extension
  !     FracRoot1stCSinkL=fraction of primary root sink strength in axis
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
  !     PSIRoot_pvr,PSIRootTurg_vr,PSIRootOSMO_vr=root total,turgor,osmotic water potential
  !     NRoot1stTipLay_raxes=layer id where primary roottip locates
  !
  !root grows into next soil layer
  IF(FGROZ.GT.0.0_r8)THEN
    RootAge_rpvr(L,NR,NZ)=1._r8
    DO NE=1,NumPlantChemElms
      RootMyco1stStrutElms_rpvr(NE,Lnext,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,Lnext,NR,NZ)+RootNetGrowthElms(NE)*FGROZ
    ENDDO

    RootProteinC_pvr(N,Lnext,NZ) = RootProteinC_pvr(N,Lnext,NZ)+&
      AMIN1(rProteinC2RootN_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmn,Lnext,NR,NZ),&
         rProteinC2RootP_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmp,Lnext,NR,NZ))
            
    PopuRootMycoC_pvr(N,Lnext,NZ)   = PopuRootMycoC_pvr(N,Lnext,NZ)+RootMyco1stStrutElms_rpvr(ielmc,Lnext,NR,NZ)
    Root1stLenPP_rpvr(Lnext,NR,NZ)  = Root1stLenPP_rpvr(Lnext,NR,NZ)+Root1stExtenzPP*FGROZ
    Root1stRadius_pvr(N,Lnext,NZ)   = Root1stRadius_pvr(N,L,NZ)
    RootCRRadius0_rpvr(Lnext,NR,NZ) = RootCRRadius0_rpvr(L,NR,NZ)
    !divide and extend into next layer
    DO NE=1,NumPlantChemElms
      XFRE(NE)                              = FracRoot1stCSinkL*RootMycoNonstElms_rpvr(NE,N,L,NZ)
      RootMycoNonstElms_rpvr(NE,N,L,NZ)     = RootMycoNonstElms_rpvr(NE,N,L,NZ)-XFRE(NE)
      RootMycoNonstElms_rpvr(NE,N,Lnext,NZ) = RootMycoNonstElms_rpvr(NE,N,Lnext,NZ)+XFRE(NE)
    ENDDO
    PSIRoot_pvr(N,Lnext,NZ)      = PSIRoot_pvr(N,L,NZ)
    PSIRootOSMO_vr(N,Lnext,NZ)   = PSIRootOSMO_vr(N,L,NZ)
    PSIRootTurg_vr(N,Lnext,NZ)   = PSIRootTurg_vr(N,L,NZ)
    NRoot1stTipLay_raxes(NR,NZ)  = MAX(NGTopRootLayer_pft(NZ),Lnext)
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine PrimeRootsExtension

!----------------------------------------------------------------------------------------------------
  subroutine PrimeRootsWithdraw(L,NR,NZ,N,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr)
  implicit none
  integer, intent(in) :: L,NR,NZ,N

  real(r8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  integer :: LL,NN,NE,NTG
  real(r8) :: XFRD,XFRW,FracRootSinkL
  real(r8) :: XFRE
  logical :: RootDepzChk
  real(r8), parameter :: scalp=0.99999_r8
  character(len=*), parameter :: subname='PrimeRootsWithdraw'

! begin_execution
  associate(                                                                &
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    CumSoilThickness_vr          => plt_site%CumSoilThickness_vr           ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZEROS2                       => plt_site%ZEROS2                        ,& !input  :threshold zero for numerical stability,[-]
    DLYR3                        => plt_site%DLYR3                         ,& !input  :vertical thickness of soil layer, [m]
    VLSoilPoreMicP_vr            => plt_soilchem%VLSoilPoreMicP_vr         ,& !input  :volume of soil layer, [m3 d-2]
    Root1stDepz_raxes            => plt_morph%Root1stDepz_raxes            ,& !input  :root layer depth, [m]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft           ,& !input  :soil layer at planting depth, [-]
    iPlantNfixType_pft           => plt_morph%iPlantNfixType_pft           ,& !input  :N2 fixation type,[-]
    Myco_pft                     => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]
    SeedDepth_pft                => plt_morph%SeedDepth_pft                ,& !input  :seeding depth, [m]
    Root1stAxesTipDepz2Surf_pft  => plt_morph%Root1stAxesTipDepz2Surf_pft  ,& !input  :plant primary depth relative to column surface, [m]         
    NumAxesPerPrimRoot_pft         => plt_morph%NumAxesPerPrimRoot_pft         ,& !output :primary root axes number, [d-2]       
    RootAge_rpvr                 => plt_morph%RootAge_rpvr                 ,& !inoput :root age,[h]
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr     ,& !inoput :root layer element primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr    => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !inoput :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr        ,& !inoput :root layer nonstructural element, [g d-2]
    RootProteinC_pvr             => plt_biom%RootProteinC_pvr              ,& !inoput :root layer protein C, [gC d-2]
    PopuRootMycoC_pvr            => plt_biom% PopuRootMycoC_pvr            ,& !inoput :root layer C, [gC d-2]
    RootNodulStrutElms_rpvr      => plt_biom%RootNodulStrutElms_rpvr       ,& !inoput :root layer nodule element, [g d-2]
    RootNodulNonstElms_rpvr      => plt_biom%RootNodulNonstElms_rpvr       ,& !inoput :root layer nonstructural element, [g d-2]
    RootGasLossDisturb_pft       => plt_bgcr%RootGasLossDisturb_pft        ,& !inoput :gaseous flux fron root disturbance, [g d-2 h-1]
    trcg_rootml_pvr              => plt_rbgc%trcg_rootml_pvr               ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr              => plt_rbgc%trcs_rootml_pvr               ,& !inoput :root aqueous content, [g d-2]
    Root1stXNumL_pvr            => plt_morph%Root1stXNumL_pvr            ,& !inoput :root layer number primary axes, [d-2]
    Root2ndXNumL_rpvr            => plt_morph%Root2ndXNumL_rpvr            ,& !inoput :root layer number axes, [d-2]
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr              ,& !inoput :root layer length primary axes, [m d-2]
    Root2ndXNum_rpvr             => plt_morph%Root2ndXNum_rpvr             ,& !inoput :root layer number secondary axes, [d-2]
    NRoot1stTipLay_raxes         => plt_morph%NRoot1stTipLay_raxes          & !output :maximum soil layer number for root axes, [-]
  )
  !     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
  !     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
  !     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
  !     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
  !
  !     NRoot1stTipLay_raxes=deepest root layer
  !     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
  !     Root1stAxesTipDepz2Surf_pft(NR,NZ)=primary root depth from soil surface
  !     CumSoilThickness_vr=depth from soil surface to layer bottom
  !     SeedDepth_pft=seeding depth
  !     FracRootSinkL=fraction of primary+secondary root sink strength in axis
  !     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
  !     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
  !     Root1stLenPP_rpvr=primary root length
  !     RootProteinC_pvr=root protein C mass
  !     WTRTD=root C mass
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root

  call PrintInfo('beg '//subname)

  D5115: DO LL=L,NGTopRootLayer_pft(NZ)+1,-1
    !test root in layer LL-1, or root axes above seeding depth
    RootDepzChk=Root1stAxesTipDepz2Surf_pft(NR,NZ).LT.CumSoilThickness_vr(LL-1) .OR. Root1stAxesTipDepz2Surf_pft(NR,NZ).LT.SeedDepth_pft(NZ)
    !EXIT if LL-1 is not active soil or has root, no withdraw

    IF(VLSoilPoreMicP_vr(LL-1)<=ZEROS2 .OR. .not. RootDepzChk)EXIT
    IF(RootSinkC_vr(N,LL).GT.ZERO4Groth_pft(NZ))THEN
      FracRootSinkL=(Root1stSink_pvr(LL,NR)+Root2ndSink_pvr(N,LL,NR))/RootSinkC_vr(N,LL)
    ELSE
      FracRootSinkL=1.0_r8
    ENDIF

    DO NE=1,NumPlantChemElms
      RootMyco1stStrutElms_rpvr(NE,LL-1,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,LL-1,NR,NZ)+RootMyco1stStrutElms_rpvr(NE,LL,NR,NZ)
      RootMyco1stStrutElms_rpvr(NE,LL,NR,NZ) = 0._r8        
    ENDDO
    Root1stLenPP_rpvr(LL-1,NR,NZ)=Root1stLenPP_rpvr(LL-1,NR,NZ)+Root1stLenPP_rpvr(LL,NR,NZ)
    Root1stLenPP_rpvr(LL,NR,NZ)  =0._r8

    D5110: DO NN=1,Myco_pft(NZ)
      
      DO NE=1,NumPlantChemElms
        RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ)+RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ)

        RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ) = 0._r8
        XFRE = FracRootSinkL*RootMycoNonstElms_rpvr(NE,NN,LL,NZ)
        XFRE = AMAX1(AMIN1(RootMycoNonstElms_rpvr(NE,NN,LL,NZ),XFRE),-RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ))*scalp

        RootMycoNonstElms_rpvr(NE,NN,LL,NZ)   = RootMycoNonstElms_rpvr(NE,NN,LL,NZ)-XFRE
        RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ) = RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ)+XFRE
      ENDDO

      XFRW = FracRootSinkL*RootProteinC_pvr(NN,L,NZ)
      XFRW = AMAX1(AMIN1(RootProteinC_pvr(NN,LL,NZ),XFRW),-RootProteinC_pvr(NN,LL-1,NZ))*scalp
      RootProteinC_pvr(NN,LL,NZ)    = RootProteinC_pvr(NN,LL,NZ)-XFRW
      RootProteinC_pvr(NN,LL-1,NZ)  = RootProteinC_pvr(NN,LL-1,NZ)+XFRW

      XFRD = FracRootSinkL*PopuRootMycoC_pvr(NN,LL,NZ)
      XFRD = AMAX1(AMIN1(PopuRootMycoC_pvr(NN,LL,NZ),XFRD),-PopuRootMycoC_pvr(NN,LL-1,NZ))*scalp
      PopuRootMycoC_pvr(NN,LL,NZ)   = PopuRootMycoC_pvr(NN,LL,NZ)-XFRD
      PopuRootMycoC_pvr(NN,LL-1,NZ) = PopuRootMycoC_pvr(NN,LL-1,NZ)+XFRD
      !
      !     WITHDRAW GASES IN PRIMARY ROOTS
      !
      !     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
      !     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
      !     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
      !     FracRootSinkL=fraction of primary root sink strength in axis
      !
      DO NTG=idg_beg,idg_NH3
        RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-FracRootSinkL &
          *(trcg_rootml_pvr(idg_CO2,NN,LL,NZ)+trcs_rootml_pvr(idg_CO2,NN,LL,NZ))
        trcg_rootml_pvr(NTG,NN,LL,NZ) = (1.0_r8-FracRootSinkL)*trcg_rootml_pvr(NTG,NN,LL,NZ)
        trcs_rootml_pvr(NTG,NN,LL,NZ) = (1.0_r8-FracRootSinkL)*trcs_rootml_pvr(NTG,NN,LL,NZ)
      ENDDO
    ENDDO D5110
    !
    !     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
    !
    !     Root2ndXNum_rpvr,Root2ndXNumL_rpvr=number of secondary root axes
    !     Root1stXNumL_pvr=number of primary root axes
    !     Root1stLenPP_rpvr=primary root length
    !     CumSoilThickness_vr=depth from soil surface to layer bottom
    !     SeedDepth_pft=seeding depth
    !
    Root2ndXNumL_rpvr(N,LL,NZ)   = AZMAX1(Root2ndXNumL_rpvr(N,LL,NZ)-Root2ndXNum_rpvr(N,LL,NR,NZ))
    Root2ndXNumL_rpvr(N,LL-1,NZ) = AZMAX1(Root2ndXNumL_rpvr(N,LL-1,NZ)+Root2ndXNum_rpvr(N,LL,NR,NZ))
    Root2ndXNum_rpvr(N,LL,NR,NZ) = 0._r8
    Root1stXNumL_pvr(LL,NZ)     = Root1stXNumL_pvr(LL,NZ)-NumAxesPerPrimRoot_pft(NZ)

    IF(LL-1.GT.NGTopRootLayer_pft(NZ))THEN
      !layer LL-1 is below the root tip layer
      Root1stLenPP_rpvr(LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness_vr(LL-1)-Root1stDepz_raxes(NR,NZ))
    ELSE
      !layer LL-1 is above root tip layer
      Root1stLenPP_rpvr(LL-1,NR,NZ)=DLYR3(LL-1)-(CumSoilThickness_vr(LL-1)-Root1stDepz_raxes(NR,NZ)) &
        -(SeedDepth_pft(NZ)-CumSoilThickness_vr(LL-2))
    ENDIF
    !
    !     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
    !
    !     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
    !     FracRootSinkL=fraction of primary root sink strength in axis
    !     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
    !     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
    !
    IF(is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
      DO NE=1,NumPlantChemElms
        XFRE = FracRootSinkL*RootNodulStrutElms_rpvr(NE,LL,NZ)
        XFRE = AMAX1(AMIN1(RootNodulStrutElms_rpvr(NE,LL,NZ),XFRE),-RootNodulStrutElms_rpvr(NE,LL-1,NZ))*scalp
        RootNodulStrutElms_rpvr(NE,LL,NZ) = RootNodulStrutElms_rpvr(NE,LL,NZ)-XFRE
        RootNodulStrutElms_rpvr(NE,LL-1,NZ) = RootNodulStrutElms_rpvr(NE,LL-1,NZ)+XFRE

        XFRE = FracRootSinkL*RootNodulNonstElms_rpvr(NE,LL,NZ)
        XFRE = AMAX1(AMIN1(RootNodulNonstElms_rpvr(NE,LL,NZ),XFRE),-RootNodulNonstElms_rpvr(NE,LL-1,NZ))*scalp
        RootNodulNonstElms_rpvr(NE,LL,NZ)   = RootNodulNonstElms_rpvr(NE,LL,NZ)-XFRE
        RootNodulNonstElms_rpvr(NE,LL-1,NZ) = RootNodulNonstElms_rpvr(NE,LL-1,NZ)+XFRE
      ENDDO
    ENDIF
    NRoot1stTipLay_raxes(NR,NZ) = MAX(NGTopRootLayer_pft(NZ),LL-1)
    RootAge_rpvr(LL-1,NR,NZ)    = 0._r8
  ENDDO D5115
  call PrintInfo('end '//subname)
  end associate
  end subroutine PrimeRootsWithdraw
!----------------------------------------------------------------------------------------------------

  function morphogen_signal(Abase, dist_scaled)result(A)
  !
  !Description:
  !compute morphogen signal at scaled distance dist_scaled
  implicit none
  real(r8), intent(in) :: Abase        !base morephogen signal, [0-1]
  real(r8), intent(in) :: dist_scaled  !scaled distance, [-]

  real(r8) :: a

  if(dist_scaled<6._r8)then
    A=Abase+(1._r8-Abase)*exp(-dist_scaled)
  else
  A=Abase
  endif
  end function morphogen_signal
!----------------------------------------------------------------------------------------------------
  subroutine SummarizeRootSink(I,J,NZ,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,flag2ndGrowth)
  !
  !DESCRIPTION
  !summarize root sink strength for growth allocation and shoot-root coupling
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8),INTENT(OUT) :: RootSinkC_vr(pltpar%jroots,JZ1)       !conductance for nonstrucal CNP allocation
  real(r8),intent(out) :: Root1stSink_pvr(JZ1,MaxNumRootAxes)
  real(r8),intent(out) :: Root2ndSink_pvr(pltpar%jroots,JZ1,MaxNumRootAxes)
  real(r8),INTENT(OUT) :: RootSinkC(pltpar%jroots)
  logical, intent(out) :: flag2ndGrowth(JZ1,MaxNumRootAxes)
  integer :: N,L,K,NR,NE,ntu,LTip
  REAL(R8) :: Root1stLocDepz_vr(NumCanopyLayers1,JZ1)
  real(r8) :: RecoRootMycoC4Nup,CUPRO,CUPRC
  real(r8) :: DistCanopyPrimeRootTip       !distance between primary root tip and canopy top, which characterizes source-sink distance
  real(r8) :: RootEffDepz,RTSKP
  real(r8) :: RTSKS,rscal,TipRadius
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: RCO2flx,dUptakeN,dUptakeP,DTransptTube,AreaTranspt,dist2Tip
  real(r8) :: DistRootEffDepz,dist_scaled,BaseSink
  real(r8), parameter :: CMassCost4NiUptk=0.86_r8  !=12/14.
  character(len=*), parameter :: subname='SummarizeRootSink'

  associate(                                                             &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft               ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantRootProfile_pft      => plt_pheno%iPlantRootProfile_pft       ,& !input  :plant growth type (vascular, non-vascular),[-]
    VLSoilPoreMicP_vr          => plt_soilchem%VLSoilPoreMicP_vr        ,& !input  :volume of soil layer, [m3 d-2]
    ZEROS2                     => plt_site%ZEROS2                       ,& !input  :threshold zero for numerical stability,[-]
    NU                         => plt_site%NU                           ,& !input  :current soil surface layer number, [-]
    CumSoilThickness_vr        => plt_site%CumSoilThickness_vr          ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZERO                       => plt_site%ZERO                         ,& !input  :threshold zero for numerical stability, [-]
    Root1stMaxRadius1_pft      => plt_morph%Root1stMaxRadius1_pft       ,& !input  :root radius primary axes, [m]        
    DLYR3                      => plt_site%DLYR3                        ,& !input  :vertical thickness of soil layer, [m]
    RootOUlmNutUptake_pvr      => plt_rbgc%RootOUlmNutUptake_pvr        ,& !input  :root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
    RootCUlmNutUptake_pvr      => plt_rbgc%RootCUlmNutUptake_pvr        ,& !input  :root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
    CanopyHeight4WatUptake_pft => plt_morph%CanopyHeight4WatUptake_pft  ,& !input  :canopy height, [m]
    RootMatureAge_pft          => plt_morph%RootMatureAge_pft           ,& !input  : Root maturation age, [h]
    RootAge_rpvr               => plt_morph%RootAge_rpvr                ,& !input  : root age,[h]
    NRoot1stTipLay_raxes       => plt_morph%NRoot1stTipLay_raxes        ,& !input  :maximum soil layer number for root axes, [-]    
    MorphogenBase_pft          => plt_pheno%MorphogenBase_pft           ,& !input  : baseline morphogen signal strength, [-]
    PSIRoot_pvr                => plt_ew%PSIRoot_pvr                    ,& !input  :root total water potential, [Mpa]    
    Myco_pft                   => plt_morph%Myco_pft                    ,& !input  :mycorrhizal type (no or yes),[-]
    Root1stRadius_rpvr         => plt_morph%Root1stRadius_rpvr          ,& !input  :root layer radius for each primary axes, [m]
    iPlantTurnoverPattern_pft  => plt_pheno%iPlantTurnoverPattern_pft   ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    Root1stDepz_raxes          => plt_morph%Root1stDepz_raxes           ,& !input  :root layer depth, [m]
    HypocotHeight_pft          => plt_morph%HypocotHeight_pft           ,& !input  :cotyledon height, [m]
    Root1stMaxRadius_pft       => plt_morph%Root1stMaxRadius_pft        ,& !input  :maximum radius of primary roots, [m]        
    Root2ndRadius_rpvr         => plt_morph%Root2ndRadius_rpvr          ,& !input  :root layer radius secondary axes, [m]
    Root2ndXNum_rpvr           => plt_morph%Root2ndXNum_rpvr            ,& !input  :root layer number secondary axes, [d-2]
    Root2ndEffLen4uptk_rpvr    => plt_morph%Root2ndEffLen4uptk_rpvr     ,& !input  :Layer effective root length four resource uptake, [m]
    SeedDepth_pft              => plt_morph%SeedDepth_pft               ,& !input  :seeding depth, [m]
    MaxSoiL4Root_pft           => plt_morph%MaxSoiL4Root_pft            ,& !input  :maximum soil layer number for all root axes,[-]
    NumPrimeRootAxes_pft       => plt_morph%NumPrimeRootAxes_pft        ,& !input  :root primary axis number,[-]
    NumAxesPerPrimRoot_pft     => plt_morph%NumAxesPerPrimRoot_pft      ,& !output :primary root axes number, [d-2]    
    Root2ndSinkWeight_pvr      => plt_morph%Root2ndSinkWeight_pvr       ,& !output :Secondary root nonst element sink profile, [d-2]
    Root1stSinkWeight_pvr      => plt_morph%Root1stSinkWeight_pvr       ,& !output :primary root nonst element sink profile, [d-2]
    RootNutUptakeN_pft         => plt_rbgc%RootNutUptakeN_pft           ,& !inoput :total N uptake by plant roots, [gN d-h2 h-1]
    RootNutUptakeP_pft         => plt_rbgc%RootNutUptakeP_pft           ,& !inoput :total P uptake by plant roots, [gP d-h2 h-1]
    RootMycoNonstElms_rpvr     => plt_biom%RootMycoNonstElms_rpvr       ,& !inoput :root layer nonstructural element, [g d-2]
    RootNutUptake_pvr          => plt_rbgc%RootNutUptake_pvr            ,& !inoput :root uptake of Nutrient band, [g d-2 h-1]
    RootCO2EmisPot_pvr         => plt_rbgc%RootCO2EmisPot_pvr           ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr         => plt_rbgc%RootRespPotent_pvr           ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootCO2Autor_pvr           => plt_rbgc%RootCO2Autor_pvr              & !inoput :root respiration constrained by O2, [g d-2 h-1]
  )

!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER
  call PrintInfo('beg '//subname)

  RootSinkC_vr    = 0._R8
  Root1stSink_pvr = 0._r8
  Root2ndSink_pvr = 0._r8
  RootSinkC       = 0._r8
  RCO2flx         = 0._r8
!   call SumRootBiome(NZ,mass_inital)
  flag2ndGrowth=.false.
  D4995: DO N=1,Myco_pft(NZ)
    Ltip=NU
    D4990: DO L=MaxSoiL4Root_pft(NZ),NU,-1
      !
      !     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
      !     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
      !
      !     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
      !     RecoRootMycoC4Nup=C respiration for nutrient uptake
      !     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
      !     why is 0.86? it refers to C cost for N asimilation

      IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
        !accumulate C cost for nitrogen (NO3,NO3B,NH4,NH4B) and phosphorus (H1PO4,H2PO4,H1PO4B,H2PO4B) uptake

        RecoRootMycoC4Nup=0._r8;CUPRO=0._r8;CUPRC=0._r8
        do ntu=ids_nutb_beg+1,ids_nuts_end
          if(ntu==ids_NO2B .or. ntu==ids_NO2)cycle
          RecoRootMycoC4Nup = RecoRootMycoC4Nup+RootNutUptake_pvr(ntu,N,L,NZ)
          CUPRO             = CUPRO+RootOUlmNutUptake_pvr(ntu,N,L,NZ)
          CUPRC             = CUPRC+RootCUlmNutUptake_pvr(ntu,N,L,NZ)
        enddo  
        !0.86 gC per gN or gP uptake 
        RecoRootMycoC4Nup = CMassCost4NiUptk*RecoRootMycoC4Nup
        CUPRO             = CMassCost4NiUptk*CUPRO
        CUPRC             = CMassCost4NiUptk*CUPRC
        !
        !
        !     RootCO2Autor_pvr=total root respiration
        !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
        !     RecoRootMycoC4Nup=C respiration for nutrient uptake
        !     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
        !     CPOOLR=non-structural C mass in root
        !
        RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+CUPRO
        RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+CUPRC

        !correct small magnitude values to avoid numerical downflow
        RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=AZERO(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))

        if(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)<RecoRootMycoC4Nup .and. RecoRootMycoC4Nup>0._r8)then
          rscal=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)/RecoRootMycoC4Nup*0.99999_r8          
          
          do ntu=ids_nutb_beg+1,ids_nuts_end
            if(ntu==ids_NO2B .or. ntu==ids_NO2)cycle
            RootNutUptake_pvr(ntu,N,L,NZ) = RootNutUptake_pvr(ntu,N,L,NZ)*rscal
          enddo
          RecoRootMycoC4Nup=RecoRootMycoC4Nup*rscal
        endif

        !root/myco non-structrual C used for root N/P uptake
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
        dUptakeN= RootNutUptake_pvr(ids_NH4,N,L,NZ) +RootNutUptake_pvr(ids_NH4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_NO3,N,L,NZ)+RootNutUptake_pvr(ids_NO3B,N,L,NZ)
        dUptakeP=RootNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootNutUptake_pvr(ids_H2PO4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ)

        RootMycoNonstElms_rpvr(ielmn,N,L,NZ)=RootMycoNonstElms_rpvr(ielmn,N,L,NZ)+dUptakeN
        RootMycoNonstElms_rpvr(ielmp,N,L,NZ)=RootMycoNonstElms_rpvr(ielmp,N,L,NZ)+dUptakeP

        RootNutUptakeN_pft(NZ)=RootNutUptakeN_pft(NZ)+dUptakeN
        RootNutUptakeP_pft(NZ)=RootNutUptakeP_pft(NZ)+dUptakeP

        D4985: DO NR=1,NumPrimeRootAxes_pft(NZ)
          !
          !     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
          !
          !     RTDP1=primary root depth from soil surface
          !     CumSoilThickness_vr=depth from soil surface to layer bottom
          !     CanopyHeight4WatUptake_pft=canopy height for water uptake
          !     RTSK=relative primary root sink strength
          !     Root1stSink_pvr=primary root sink strength
          !     NumAxesPerPrimRoot_pft=number of primary root axes
          !     RRAD1,Root2ndRadius_rpvr=primary, secondary root radius
          !     RootSinkC,RootSinkC_vr=total root sink strength
          !
          IF(N.EQ.ipltroot)THEN
            Root1stLocDepz_vr(NR,L) = AZMAX1(Root1stDepz_raxes(NR,NZ)-CumSoilThickness_vr(L-1)-RTDPX)
            Root1stLocDepz_vr(NR,L) = AZMAX1(AMIN1(DLYR3(L),Root1stLocDepz_vr(NR,L))-AZMAX1(SeedDepth_pft(NZ)-CumSoilThickness_vr(L-1)-HypocotHeight_pft(NZ)))
            RootEffDepz             = AMAX1(SeedDepth_pft(NZ),CumSoilThickness_vr(L-1))+0.5_r8*Root1stLocDepz_vr(NR,L)
            DistRootEffDepz         = RootEffDepz+CanopyHeight4WatUptake_pft(NZ)

            DistCanopyPrimeRootTip = Root1stDepz_raxes(NR,NZ)+CanopyHeight4WatUptake_pft(NZ)
            DTransptTube           = AMIN1(ZSTX,AMAX1(FSTK*Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(N,NZ)))
            AreaTranspt            = 2._r8*AMAX1(Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(N,NZ))*DTransptTube-DTransptTube**2

            !pi*(R^2-(R-d)^2)=pi*(2*R*d-d^2)

            IF(Root1stDepz_raxes(NR,NZ).GT.CumSoilThickness_vr(L-1))THEN              
              IF(Root1stDepz_raxes(NR,NZ).LT.CumSoilThickness_vr(L))THEN
                !Root tip in layer L
                Ltip=L
                TipRadius             = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
                Root1stSink_pvr(L,NR) = RTSK(iPlantRootProfile_pft(NZ))*NumAxesPerPrimRoot_pft(NZ)*TipRadius**2/DistCanopyPrimeRootTip
                RootSinkC_vr(N,L)      = RootSinkC_vr(N,L)+Root1stSink_pvr(L,NR)
              elseif(iPlantTurnoverPattern_pft(NZ).NE.0 .and. is_plant_treelike(iPlantRootProfile_pft(NZ)) .and. &
                RootAge_rpvr(L,NR,NZ)>RootMatureAge_pft(NZ) .and. RootEffDepz<Root1stDepz_raxes(NR,NZ) .and. &
                Root1stDepz_raxes(NR,NZ)>7.5_r8*TipRadius)THEN
                !non-tip layer
                flag2ndGrowth(L,NR)=.true.
                if(Root1stSink_pvr(MaxSoiL4Root_pft(NZ),NR).LE.ZERO .and. NRoot1stTipLay_raxes(NR,NZ)==MaxSoiL4Root_pft(NZ))then
                  !roottip stop going deeper at the bottom, and 
                  Ltip      = MaxSoiL4Root_pft(NZ)
                  TipRadius = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,Ltip,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
                  BaseSink  = RTSK(iPlantRootProfile_pft(NZ))*NumAxesPerPrimRoot_pft(NZ)*TipRadius**2/DistCanopyPrimeRootTip
                else
                  TipRadius = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,Ltip,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
                  BaseSink  = Root1stSink_pvr(Ltip,NR)
                endif
                dist_scaled           = (Root1stDepz_raxes(NR,NZ)-RootEffDepz)/(2._r8*TipRadius)
                Root1stSink_pvr(L,NR) = BaseSink*morphogen_signal(MorphogenBase_pft(NZ), dist_scaled)*AreaTranspt/TipRadius**2
                RootSinkC_vr(N,L)     = RootSinkC_vr(N,L)+Root1stSink_pvr(L,NR)
              ENDIF
            ENDIF

            !
            !     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
            !     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
            !
            !     Root1stLocDepz_vr=depth of primary root axis in layer L, local bottom is DLYR3(L)
            !     RTDP1=primary root depth from soil surface
            !     CumSoilThickness_vr=depth from soil surface to layer bottom
            !     RTDPX=distance behind growing point for secondary roots
            !     DLYR=layer thickness
            !     SeedDepth_pft=seeding depth
            !     HypocotHeight_pft=hypocotyledon height is the vertical distance from the base of a seedling’s 
            !      stem—just above where the root (radicle) ends—to the point where the cotyledons (seed leaves) are attached
            !     CanopyHeight4WatUptake_pft=canopy height for water uptake
            !     RootEffDepz=secondary root depth from canopy
            !     RTSKP,RTSKS=primary,secondary root sink strength
            !     RTN2=number of secondary root axes
            !     Root2ndSink_pvr=total secondary root sink strength
            !     Root2ndEffLen4uptk_rpvr=average secondary root length
            !     RootSinkC,RootSinkC_vr=total root sink strength
            !

            IF(DistRootEffDepz.GT.ZERO)THEN
              !In the Münch model, it is assumed the actual phloem flow is via a collection of thin sieve pores (about 1 um radius), and the R**2 rule counts
              !number of pores. 
              RTSKP = NumAxesPerPrimRoot_pft(NZ)*AreaTranspt/DistRootEffDepz                
              RTSKS = safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)*Root2ndRadius_rpvr(N,L,NZ)**2,Root2ndEffLen4uptk_rpvr(N,L,NZ))
               
              IF(RTSKP+RTSKS.GT.ZERO4Groth_pft(NZ))THEN
                Root2ndSink_pvr(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
              ELSE
                Root2ndSink_pvr(N,L,NR)=0._r8
              ENDIF
            ELSE
              Root2ndSink_pvr(N,L,NR)=0._r8
            ENDIF
          ELSE
            !mycorrhizae
            Root2ndSink_pvr(N,L,NR)=safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)*Root2ndRadius_rpvr(N,L,NZ)**2,Root2ndEffLen4uptk_rpvr(N,L,NZ))
          ENDIF

          RootSinkC_vr(N,L) = RootSinkC_vr(N,L)+Root2ndSink_pvr(N,L,NR)
        ENDDO D4985        
        RootSinkC(N)      = RootSinkC(N)+RootSinkC_vr(N,L)
      ENDIF
    ENDDO D4990
    Root2ndSinkWeight_pvr(:,N,NZ)=0._r8
    if(RootSinkC(N)>0._r8)then
      DO L=NU,MaxSoiL4Root_pft(NZ)
        DO NR=1,NumPrimeRootAxes_pft(NZ)
          Root2ndSinkWeight_pvr(L,N,NZ)=Root2ndSinkWeight_pvr(L,N,NZ)+Root2ndSink_pvr(N,L,NR)
        ENDDO
        Root2ndSinkWeight_pvr(L,N,NZ)=Root2ndSinkWeight_pvr(L,N,NZ)/RootSinkC(N)
      ENDDO
    endif
  ENDDO D4995

  Root1stSinkWeight_pvr(:,NZ)=0._r8
  if(RootSinkC(ipltroot)>0._r8)then
    DO L=NU,MaxSoiL4Root_pft(NZ)
      DO NR=1,NumPrimeRootAxes_pft(NZ)
        Root1stSinkWeight_pvr(L,NZ)=Root1stSinkWeight_pvr(L,NZ)+Root1stSink_pvr(L,NR)
      ENDDO
      Root1stSinkWeight_pvr(L,NZ)=Root1stSinkWeight_pvr(L,NZ)/RootSinkC(ipltroot)
    ENDDO   
  endif   
!   call SumRootBiome(NZ,mass_finale)
  call PrintInfo('end '//subname)
  end associate
  end subroutine SummarizeRootSink

!----------------------------------------------------------------------------------------------------
  subroutine RootCheck(I,J,NZ,info)
  !
  !root axes are sorted by root order and then also by vertical layers
  !check primary root, with and without vertical distribution
  implicit none
  integer , intent(in) :: I,J,NZ
  character(len=*), intent(in) :: info
  integer :: NR,L,N,NE
  real(r8) :: dRootMyco1stElm_raxs(NumPlantChemElms)

  associate(                                                          &
    NumPrimeRootAxes_pft      => plt_morph%NumPrimeRootAxes_pft      ,& !input  :root primary axis number,[-]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    NL                        => plt_site%NL                         ,& !input  :lowest soil layer number,[-]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !input  :root layer element primary axes, [g d-2]
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs         & !input  :root C primary axes, [g d-2]
  )
  DO NR=1,NumPrimeRootAxes_pft(NZ)
    
    dRootMyco1stElm_raxs(:)=0._r8
    DO L=NU,NL
      DO NE=1,NumPlantChemElms
        dRootMyco1stElm_raxs(NE)=dRootMyco1stElm_raxs(NE)+RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)  
      ENDDO
    ENDDO
    write(111,*)trim(info),' rootcheck:',I+J/24.,NZ,N,NR,&
        dRootMyco1stElm_raxs(:),RootMyco1stElm_raxs(:,NR,NZ),dRootMyco1stElm_raxs(:)-RootMyco1stElm_raxs(:,NR,NZ)
    dRootMyco1stElm_raxs(:)=dRootMyco1stElm_raxs(:)-RootMyco1stElm_raxs(:,NR,NZ)   
    if(any(abs(dRootMyco1stElm_raxs(:))>1.e-8_r8))stop
    
  ENDDO
  end associate
  end subroutine RootCheck
!----------------------------------------------------------------------------------------------------

  subroutine UpdateRootElongationZoneAge(I,J,NR,NZ,ElongateRate)
  !
  !Description
  !
  implicit none
  integer, intent(in) :: I,J,NR,NZ
  real(r8), intent(in)    :: ElongateRate  !Root length Growth rate (>0) or Retreat rate (<0)
  real(r8), parameter :: DLength=0.05_r8  !segment length [m]  
  real(r8), parameter :: dt=1._r8
  character(len=*), parameter :: subname='UpdateRootElongationZoneAge'
  real(r8) :: growth, retreat, current_len, current_age
  real(r8) :: new_total_len, new_mean_age, excess_len
  integer  :: k, curr_idx
  
  associate(                                                     &
    RootSegAges_raxes      => plt_morph%RootSegAges_raxes      , & !inoput: age of different active root segments
    NActiveRootSegs_raxes  => plt_morph%NActiveRootSegs_raxes  , & !inoput: number of active root segments
    RootSegBaseDepth_raxes => plt_morph%RootSegBaseDepth_raxes , & !inoput: base depth of different root axes
    NMaxRootSegs           => pltpar%NMaxRootSegs              , & !input:  Maximum number of root segments for age tracking
    RootSeglengths_raxes   => plt_morph%RootSeglengths_raxes   , & !inoput: root length in each segment
    IndRootSegBase_raxes   => plt_morph%IndRootSegBase_raxes   , & !inoput: Index for base segment under tracking
    IndRootSegTip_raxes    => plt_morph%IndRootSegTip_raxes      & !inoput: Index for tip segment under tracking
  )

  call PrintInfo('beg '//subname)
  ! ---------------------------------------------------------
  ! 1. AGING STEP
  ! ---------------------------------------------------------
  if (NActiveRootSegs_raxes(NR,NZ) > 0) then
    curr_idx = IndRootSegBase_raxes(NR,NZ)
    do k = 1, NActiveRootSegs_raxes(NR,NZ)
      RootSegAges_raxes(curr_idx,NR,NZ) = RootSegAges_raxes(curr_idx,NR,NZ) + dt
      
      ! Circular Step (1-based logic):
      ! ((Current - 1 + 1) % NMaxRootSegs) + 1
      curr_idx = mod(curr_idx, NMaxRootSegs) + 1
    end do
  end if

  ! ---------------------------------------------------------
  ! 2. HANDLE GROWTH (Positive Rate)
  ! ---------------------------------------------------------
  if (ElongateRate > 0.0_r8) then
    growth = ElongateRate * dt
    
    ! Initialize if empty
    if (NActiveRootSegs_raxes(NR,NZ) == 0) then
        NActiveRootSegs_raxes(NR,NZ)  = 1
        IndRootSegBase_raxes(NR,NZ)   = 1
        IndRootSegTip_raxes(NR,NZ)    = 1
        RootSeglengths_raxes(1,NR,NZ) = 0.0_r8
        RootSegAges_raxes(1,NR,NZ)    = 0.0_r8
    end if
      
    ! Get current tip data
    current_len = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)
    current_age = RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)
    
    ! Weighted Mean Age Logic
    new_total_len = current_len + growth
    new_mean_age  = (current_age * current_len) / new_total_len
    
    ! Update Tip
    RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = new_total_len
    RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)    = new_mean_age
      
    ! Check for Overflow (Segment > DLength)
    if (RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) > DLength) then
      excess_len = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) - DLength
      
      ! Cap the current full segment
      RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = DLength
      
      ! Move Tip Forward (1-based circular)
      IndRootSegTip_raxes(NR,NZ) = mod(IndRootSegTip_raxes(NR,NZ), NMaxRootSegs) + 1
      
      ! Check for Buffer Full
      if (NActiveRootSegs_raxes(NR,NZ) < NMaxRootSegs) then
        NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ) + 1
      else
        ! Buffer Full: Tip overwrites Base
        ! Base must move forward to escape
        ! Update RootSegBaseDepth_raxes(NR) because we lose the geometry of the old base
        RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ) + RootSeglengths_raxes(IndRootSegBase_raxes(NR,NZ),NR,NZ)
        
        IndRootSegBase_raxes(NR,NZ) = mod(IndRootSegBase_raxes(NR,NZ), NMaxRootSegs) + 1
      end if
      
      ! Initialize the new tip
      RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = excess_len
      RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)    = 0.0_r8
    end if
      
  ! ---------------------------------------------------------
  ! 3. HANDLE RETREAT (Negative Rate)
  ! ---------------------------------------------------------
  else if (ElongateRate < 0.0_r8) then
    retreat = abs(ElongateRate) * dt
    
    do while (retreat > 0.0_r8 .and. NActiveRootSegs_raxes(NR,NZ) > 0)
      current_len = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)
      
      if (retreat >= current_len) then
        ! Consume entire segment
        retreat = retreat - current_len
        RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = 0.0_r8
        RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)    = 0.0_r8
        
        if (NActiveRootSegs_raxes(NR,NZ) > 1) then
          ! Move Tip Backward (1-based circular)
          ! Logic: ((Current - 1 - 1 + N) % N) + 1
          ! Simplified: ((Current - 2 + N) % N) + 1
          IndRootSegTip_raxes(NR,NZ) = mod(IndRootSegTip_raxes(NR,NZ) - 2 + NMaxRootSegs, NMaxRootSegs) + 1
          NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ) - 1
        else
          ! Root is gone
          NActiveRootSegs_raxes(NR,NZ) = 0
          retreat                   = 0.0_r8
        end if
      else
        ! Partial consumption
        RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) - retreat
        retreat = 0.0_r8
      end if
    end do
  end if
  call PrintInfo('end '//subname)
  end associate      
  end subroutine UpdateRootElongationZoneAge
!----------------------------------------------------------------------------------------------------

  subroutine FindMatureRootSegs(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  character(len=*), parameter :: subname='FindMatureRootSegs'
  real(r8), parameter :: DLength=0.05_r8  !segment length [m]    
  real(r8), parameter :: RootLminMat=0.05_r8 !minimum length between growing tip and zone of secondary growth [m]
  INTEGER :: k,NR,KD
  logical :: foundMatureSeg

  associate(                                                      &
    NumPrimeRootAxes_pft    => plt_morph%NumPrimeRootAxes_pft   , & !input: root primary axis number,[-]  
    RootSegAges_raxes       => plt_morph%RootSegAges_raxes      , & !input: age of different active root segments
    NActiveRootSegs_raxes   => plt_morph%NActiveRootSegs_raxes  , & !inoput: number of active root segments
    RootSegBaseDepth_raxes  => plt_morph%RootSegBaseDepth_raxes , & !inoput: base depth of different root axes, [m]
    Root1stDepz_raxes       => plt_morph%Root1stDepz_raxes      , & !input : root layer depth, [m]    
    NMaxRootSegs            => pltpar%NMaxRootSegs              , & !input : Maximum number of root segments for age tracking
    RootMatureAge_pft       => plt_morph%RootMatureAge_pft      , & !input : Root maturation age, [h]
    SeedDepth_pft           => plt_morph%SeedDepth_pft          , & !input  :seeding depth, [m]
    RootSeglengths_raxes    => plt_morph%RootSeglengths_raxes   , & !inoput: root length in each segment, [m]
    IndRootSegBase_raxes    => plt_morph%IndRootSegBase_raxes   , & !inoput: Index for base segment under tracking
    IndRootSegTip_raxes     => plt_morph%IndRootSegTip_raxes      & !inoput: Index for tip segment under tracking
  )
  call PrintInfo('beg '//subname)
  DO NR=1,NumPrimeRootAxes_pft(NZ)
    foundMatureSeg=.false.

    !the tracking array is of shape [1,2,...->tip,...,base->,...,max]
    IF(IndRootSegBase_raxes(NR,NZ) > IndRootSegTip_raxes(NR,NZ))then
      DO K=IndRootSegTip_raxes(NR,NZ),1,-1
        if(RootSegAges_raxes(K,NR,NZ)>=RootMatureAge_pft(NZ))then
          foundMatureSeg=.true.
          KD=K-1
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+DLength*KD
          RootSegBaseDepth_raxes(NR,NZ) = min(RootSegBaseDepth_raxes(NR,NZ)+RootSeglengths_raxes(K,NR,NZ),Root1stDepz_raxes(NR,NZ))
          IF(Root1stDepz_raxes(NR,NZ)<RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            

          RootSeglengths_raxes(1:K,NR,NZ) = 0._r8
          RootSegAges_raxes(1:K,NR,NZ)    = 0._r8
          IndRootSegBase_raxes(NR,NZ)     = K+1
          NActiveRootSegs_raxes(NR,NZ)    = NActiveRootSegs_raxes(NR,NZ)-KD-1
          exit
        ENDIF
      ENDDO
      if(foundMatureSeg)then
        if(IndRootSegBase_raxes(NR,NZ)>IndRootSegTip_raxes(NR,NZ))then
          !the elongation zone is dead.
          NActiveRootSegs_raxes(NR,NZ)=0
          RootSegBaseDepth_raxes(NR,NZ)=min(RootSegBaseDepth_raxes(NR,NZ),Root1stDepz_raxes(NR,NZ))
          IF(Root1stDepz_raxes(NR,NZ)<RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            
        endif
        cycle      
      endif  

      DO K=NMaxRootSegs,IndRootSegBase_raxes(NR,NZ),-1
        if(RootSegAges_raxes(K,NR,NZ)>=RootMatureAge_pft(NZ))then
          foundMatureSeg                = .true.
          Kd=(K-IndRootSegBase_raxes(NR,NZ))
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+DLength*Kd
          RootSegBaseDepth_raxes(NR,NZ) = min(RootSegBaseDepth_raxes(NR,NZ)+RootSeglengths_raxes(K,NR,NZ),Root1stDepz_raxes(NR,NZ))
          IF(Root1stDepz_raxes(NR,NZ)<RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            

          RootSegAges_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ)    = 0._r8
          RootSeglengths_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ) = 0._r8

          IndRootSegBase_raxes(NR,NZ)  = K+1
          NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ)-KD-1

          exit  
        endif
      ENDDO
    else
      DO K=IndRootSegTip_raxes(NR,NZ),IndRootSegBase_raxes(NR,NZ),-1
        if(RootSegAges_raxes(K,NR,NZ)>=RootMatureAge_pft(NZ))then
          foundMatureSeg=.true.
          KD=K-IndRootSegBase_raxes(NR,NZ)
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+DLength*KD          
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+RootSeglengths_raxes(K,NR,NZ)
          IF(Root1stDepz_raxes(NR,NZ)<RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            

          RootSegAges_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ)    = 0._r8
          RootSeglengths_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ) = 0._r8

          NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ)-KD-1
          IndRootSegBase_raxes(NR,NZ)  = K+1
          exit
        ENDIF
      ENDDO
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine FindMatureRootSegs  
!----------------------------------------------------------------------------------------------------
  subroutine ComputeSoilStress(I,J,Pstress_vr)
  implicit none
  INTEGER, intent(in) :: I,J
  real(r8), intent(out) :: Pstress_vr(JZ1)     !Effective Vertical Stress
  integer :: L
  character(len=*), parameter :: subname='ComputeSoilStress'
  associate(                                                     &
    NU                         => plt_site%NU                  , & !input : current soil surface layer number, [-]
    SoilWeightStress_vr        => plt_site%SoilWeightStress_vr , & !input : soil weight stress on root thickening, [MPa]
    SoilSuctStress_vr          => plt_site%SoilSuctStress_vr   , & !input : soil suction stress on root thickening, [MPa]
    rSat_vr                    => plt_site%rSat_vr               & !input : relative soil saturation, [-]            
  )
  call PrintInfo('beg '//subname)

  Pstress_vr(1)=SoilWeightStress_vr(1)*0.5_r8+rSat_vr(1)*abs(SoilSuctStress_vr(1))

  DO L=2,JZ1
    Pstress_vr(L)=Pstress_vr(L-1)+(SoilWeightStress_vr(L-1)+SoilWeightStress_vr(L))*0.5_r8+rSat_vr(L)*abs(SoilSuctStress_vr(L))
  enddo
  call PrintInfo('end '//subname)

  end associate
  end subroutine ComputeSoilStress
!----------------------------------------------------------------------------------------------------

  ![tail]
end module RootMod
