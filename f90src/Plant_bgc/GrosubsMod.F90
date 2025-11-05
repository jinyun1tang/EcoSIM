module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod,         only: safe_adb, AZMAX1,real_truncate
  use data_kind_mod,       only: r8 => DAT_KIND_R8
  use EcoSIMCtrlMod,       only: lverb
  use EcoSiMParDataMod,    only: pltpar
  use RootMod,             only: RootBGCModel
  use PlantNonstElmDynMod, only: PlantNonstElmTransfer
  use PlantDebugMod,       only: PrintRootTracer
  use DebugToolMod,        only: PrintInfo
  use EcosimConst
  use PlantBGCPars
  use PlantAPIData
  use PhotoSynsMod
  use PlantMathFuncMod
  use LitterFallMod
  use PlantBranchMod
  use PlantBalMod
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
! DIMENSION VCO2(400,366,05)
!
!     RTSK=relative primary root sink strength 0.25=shallow,4.0=deep root profile
!     FXRN=rate constant for plant-bacteria nonstructl C,N,P exchange (h-1)
!     RateK4ShootSeaStoreNonstEXfer=rate constant for leaf-storage nonstructl C,N,P exchange (h-1)
!     RateK4RootSeaStorNonstEXfer=rate constant for root-storage nonstructl C,N,P exchange (h-1)
!     FPART=parameter for allocating nonstructural C to shoot organs
!     FXSH,FXRT=shoot-root partitioning of storage C during leafout
!     FRSV=rate constant for remobiln of storage C,N,P during leafout (h-1)
!     FXFY,FXFZ=rate constant for leaf-reserve nonstructural C,N,P exchange (h-1)
!     EFIRE=combustion  of N,P relative to C
!     PSILY=canopy water potential below which leafoff is induced (MPa)
!     HourFailGrainFill_brchY=number of hours after physiol maturity required for senescence
!     HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
!     GVMX=specific oxidation rate of nonstructural C during leafout at 25 C
!

  integer  :: curday,curhour
  public :: GrowPlants

  contains
  ![header]

!----------------------------------------------------------------------------------------------------
  subroutine GrowPlants(I,J)
!
!     THIS subroutine CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
!
  use PlantDisturbsMod, only : RemoveBiomassByDisturbance,StageDisturbances
  implicit none
  integer, intent(in) :: I, J

  real(r8) :: CanopyHeight_copy(JP1)
  integer :: L,K,M
  integer :: NZ,NE
  real(r8) :: tvegE(NumPlantChemElms),tvegE1(NumPlantChemElms)
! begin_execution
  associate(                                                               &
    IsPlantActive_pft            => plt_pheno%IsPlantActive_pft           ,& !input  :flag for living pft, [-]
    NP                           => plt_site%NP                           ,& !input  :current number of plant species,[-]
    NP0                          => plt_site%NP0                          ,& !input  :intitial number of plant species,[-]
    CanopyHeight_pft             => plt_morph%CanopyHeight_pft             & !inoput :canopy height, [m]
  )
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
!
!     INITIALIZE SENESCENCE ARRAYS
  DO NZ=1,NP0
    CanopyHeight_copy(NZ)                  = CanopyHeight_pft(NZ)
    CanopyHeight_pft(NZ)                   = 0._r8
    plt_rbgc%canopy_growth_pft(NZ)         = 0._r8
    plt_biom%RootMycoMassElm_pvr(:,:,:,NZ) = 0._r8

  ENDDO  
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
  D9985: DO NZ=1,NP
    call StageDisturbances(I,J,NZ)

    IF(IsPlantActive_pft(NZ).EQ.iActive .and. plt_site%PlantPopulation_pft(NZ)>plt_site%ZEROS)THEN      

      call GrowOnePlant(I,J,NZ,CanopyHeight_copy)

    ENDIF  

    call AccumulateStates(I,J,NZ)

    call RemoveBiomassByDisturbance(I,J,NZ)

  ENDDO D9985
!
  call LiveDeadTransformation(I,J)
  
  DO NZ=1,NP

    call SumPlantBiome(I,J,NZ,'exgrosubs')

  ENDDO

  end associate
  END subroutine GrowPlants

!----------------------------------------------------------------------------------------------------
  subroutine LiveDeadTransformation(I,J)
  implicit none
  integer, intent(in) :: I,J

  integer :: L,K,NZ,M,NE,NB
  real(r8) :: XFRC,XFRN,XFRP,XFRE
  real(r8), parameter :: StandingDeadKd=1.5814E-05_r8  !first-order decay rate of standing dead biomass to litter
!     begin_execution

  associate(                                                                    &
    CH4ByFire_CumYr_pft            => plt_distb%CH4ByFire_CumYr_pft            ,& !input  :plant CH4 emission from fire, [g d-2 ]
    CO2ByFire_CumYr_pft            => plt_distb%CO2ByFire_CumYr_pft            ,& !input  :plant CO2 emission from fire, [g d-2 ]
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr             ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    EcoHavstElmntCum_pft           => plt_distb%EcoHavstElmntCum_pft           ,& !input  :total plant element harvest, [gC d-2 ]
    GrossCO2Fix_pft                => plt_bgcr%GrossCO2Fix_pft                 ,& !input  :total gross CO2 fixation, [gC d-2 ]
    GrossResp_pft                  => plt_bgcr%GrossResp_pft                   ,& !input  :total plant respiration, [gC d-2 ]
    HourReq4LeafOut_brch           => plt_pheno%HourReq4LeafOut_brch           ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    Hours4Leafout_brch             => plt_pheno%Hours4Leafout_brch             ,& !input  :heat requirement for spring leafout/dehardening, [h]
    IsPlantActive_pft              => plt_pheno%IsPlantActive_pft              ,& !input  :flag for living pft, [-]
    MaxNumRootLays                 => plt_site%MaxNumRootLays                  ,& !input  :maximum root layer number,[-]
    N2ObyFire_CumYr_pft            => plt_distb%N2ObyFire_CumYr_pft            ,& !input  :plant N2O emission from fire, [g d-2 ]
    NH3Emis_CumYr_pft              => plt_bgcr%NH3Emis_CumYr_pft               ,& !input  :total canopy NH3 flux, [gN d-2 ]
    NH3byFire_CumYr_pft            => plt_distb%NH3byFire_CumYr_pft            ,& !input  :plant NH3 emission from fire, [g d-2 ]
    NP0                            => plt_site%NP0                             ,& !input  :intitial number of plant species,[-]
    NetCumElmntFlx2Plant_pft       => plt_pheno%NetCumElmntFlx2Plant_pft       ,& !input  :effect of canopy element status on seed set, [-]
    NumOfBranches_pft              => plt_morph%NumOfBranches_pft              ,& !input  :number of branches,[-]
    PO4byFire_CumYr_pft            => plt_distb%PO4byFire_CumYr_pft            ,& !input  :plant PO4 emission from fire, [g d-2 ]
    PlantExudElm_CumYr_pft         => plt_rbgc%PlantExudElm_CumYr_pft          ,& !input  :total net root element uptake (+ve) - exudation (-ve), [gC d-2 ]
    PlantN2Fix_CumYr_pft           => plt_bgcr%PlantN2Fix_CumYr_pft            ,& !input  :total plant N2 fixation, [g d-2 ]
    RootElms_pft                   => plt_biom%RootElms_pft                    ,& !input  :plant root element mass, [g d-2]
    SeasonalNonstElms_pft          => plt_biom%SeasonalNonstElms_pft           ,& !input  :plant stored nonstructural element at current step, [g d-2]
    ShootElms_pft                  => plt_biom%ShootElms_pft                   ,& !input  :canopy shoot structural chemical element mass, [g d-2]
    doPlantLeafOut_brch            => plt_pheno%doPlantLeafOut_brch            ,& !input  :branch phenology flag, [-]
    fTCanopyGroth_pft              => plt_pheno%fTCanopyGroth_pft              ,& !input  :canopy temperature growth function, [-]
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft          ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft      ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iYearCurrent                   => plt_site%iYearCurrent                    ,& !input  :current year,[-]
    k_fine_litr                    => pltpar%k_fine_litr                       ,& !input  :fine litter complex id
    k_woody_litr                   => pltpar%k_woody_litr                      ,& !input  :woody litter complex id
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft      ,& !inoput :total plant element LitrFall, [g d-2 ]
    LitrfallElms_pft               => plt_bgcr%LitrfallElms_pft                ,& !inoput :plant element LitrFall, [g d-2 h-1]
    LitrfallElms_pvr               => plt_bgcr%LitrfallElms_pvr                ,& !inoput :plant LitrFall element, [g d-2 h-1]
    StandDeadKCompElms_pft         => plt_biom%StandDeadKCompElms_pft          ,& !inoput :standing dead element fraction, [g d-2]
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft  ,& !inoput :total surface LitrFall element, [g d-2]
    doInitPlant_pft                => plt_pheno%doInitPlant_pft                ,& !inoput :PFT initialization flag:0=no,1=yes,[-]
    NetPrimProduct_pft             => plt_bgcr%NetPrimProduct_pft              ,& !output :total net primary productivity, [gC d-2]
    PlantinDepz_pft                => plt_morph%PlantinDepz_pft                ,& !output :planting depth, [m]
    iDayPlanting_pft               => plt_distb%iDayPlanting_pft               ,& !output :day of planting,[-]
    iYearPlanting_pft              => plt_distb%iYearPlanting_pft               & !output :year of planting,[-]
  )
  D9975: DO NZ=1,NP0
    !
    !     ACTIVATE DORMANT SEEDS
    !
    D205: DO NB=1,NumOfBranches_pft(NZ)
      IF(doInitPlant_pft(NZ).EQ.itrue)THEN
        IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
          iDayPlanting_pft(NZ)  = I
          iYearPlanting_pft(NZ) = iYearCurrent
          PlantinDepz_pft(NZ)   = 0.005_r8+CumSoilThickness_vr(0)
          doInitPlant_pft(NZ)   = ifalse   !mark plant as initialized
          exit  
        ENDIF
      ENDIF
    ENDDO D205
    !
    !     LitrFall FROM STANDING DEAD (i.e. loss of standing dead)
    !
    !     XFRC,XFRN,XFRP=LitrFall from standing dead
    !     fTCanopyGroth_pft=temperature function for canopy growth
    !    
    D6235: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        XFRE=StandingDeadKd*fTCanopyGroth_pft(NZ)*StandDeadKCompElms_pft(NE,M,NZ)
        
        IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. .not.is_plant_treelike(iPlantRootProfile_pft(NZ)))THEN
          !grass/bryophyte
          LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_litr,0,NZ)+XFRE
        ELSE
          !trees
          LitrfallElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_litr,0,NZ)+XFRE
        ENDIF
        StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ)-XFRE
      ENDDO
    ENDDO D6235
!
!     ACCUMULATE TOTAL SURFACE, SUBSURFACE LitrFall
!
!   diagnose surface literfall 
    DO K=1,pltpar%NumOfPlantLitrCmplxs      
      D6431: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          SurfLitrfalStrutElms_CumYr_pft(NE,NZ)=SurfLitrfalStrutElms_CumYr_pft(NE,NZ)+LitrfallElms_pvr(NE,M,K,0,NZ)
        ENDDO
      ENDDO D6431  
    ENDDO

!   the following purposely starts from layer 0.
    LitrfallElms_pft(1:NumPlantChemElms,NZ)=0._r8
    D8955: DO L=0,MaxNumRootLays
      DO K=1,pltpar%NumOfPlantLitrCmplxs      
        D6430: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfallElms_pft(NE,NZ)=LitrfallElms_pft(NE,NZ)+LitrfallElms_pvr(NE,M,K,L,NZ)
          enddo         
        ENDDO D6430      
      ENDDO
    ENDDO D8955    

    DO NE=1,NumPlantChemElms
      LitrfalStrutElms_CumYr_pft(NE,NZ)=LitrfalStrutElms_CumYr_pft(NE,NZ)+LitrfallElms_pft(NE,NZ)
    ENDDO

    !   GrossResp_pft < 0 respired into atmosphere
    NetPrimProduct_pft(NZ) = GrossCO2Fix_pft(NZ)+GrossResp_pft(NZ)

  ENDDO D9975
  end associate
  end subroutine LiveDeadTransformation

!----------------------------------------------------------------------------------------------------
  subroutine GrowOnePlant(I,J,NZ,CanopyHeight_copy)
  !
  !Description
  !plant growth
  use PlantDisturbsMod, only : RemoveBiomByMgmt
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  character(len=*), parameter :: subname='GrowOnePlant'
  real(r8)  :: CanopyN2Fix_pft(JP1)
  integer  :: NB,NE,KK
  integer  :: BegRemoblize
  real(r8) :: TFN6_vr(JZ1)
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT    !main branch growth allocated to leaf and petiole
  real(r8) :: RootPrimeAxsNum
  real(r8) :: TFN5
  real(r8) :: WaterStress4Groth
  real(r8) :: Stomata_Stress
  real(r8) :: TurgEff4LeafPetolExpansion,TurgEff4CanopyResp
  real(r8) :: RootSinkC_vr(pltpar%jroots,JZ1),RootSinkC(pltpar%jroots)
  real(r8) :: Root1stSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: tvegE_beg(NumPlantChemElms)
  real(r8) :: tvegE_end(NumPlantChemElms)
  real(r8) :: err,arr(12),arr1(11),PTRTM
! begin_execution
  associate(                                                 &
    NumOfBranches_pft    => plt_morph%NumOfBranches_pft     ,& !input  :number of branches,[-]
    iPlantRootState_pft  => plt_pheno%iPlantRootState_pft   ,& !input  :flag to detect root system death,[-]
    MainBranchNum_pft    => plt_morph%MainBranchNum_pft     ,& !input  :number of main branch,[-]    
    iPlantShootState_pft => plt_pheno%iPlantShootState_pft   & !input  :flag to detect canopy death,[-]
  )

  call PrintInfo('beg '//subname)
  
  IF(iPlantShootState_pft(NZ).EQ.iLive .OR. iPlantRootState_pft(NZ).EQ.iLive)THEN
    BegRemoblize        = 0
    
    call StagePlantForGrowth(I,J,NZ,TFN6_vr,CNLFW,CPLFW,&
      CNSHW,CPSHW,CNRTW,CPRTW,RootPrimeAxsNum,TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4LeafPetolExpansion,TurgEff4CanopyResp)
!
!     CALCULATE GROWTH OF EACH BRANCH


    DO  NB=1,NumOfBranches_pft(NZ)
      call GrowOneBranch(I,J,NB,NZ,TFN6_vr,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4LeafPetolExpansion,TurgEff4CanopyResp,PTRTM,BegRemoblize)
      IF(NB.EQ.MainBranchNum_pft(NZ))PTRT=PTRTM
    ENDDO

    call RootBGCModel(I,J,NZ,BegRemoblize,TFN6_vr,CNRTW,CPRTW,RootPrimeAxsNum, &
      RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)

    call PlantNonstElmTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,BegRemoblize)

    call ComputeTotalBiom(I,J,NZ)
  ENDIF

  call RemoveBiomByMgmt(I,J,NZ)  
  !
  !   RESET DEAD BRANCHES
  call ResetDeadPlant(I,J,NZ)
  
  call PrintInfo('end '//subname)  
  end associate
  end subroutine GrowOnePlant

!----------------------------------------------------------------------------------------------------
  subroutine StagePlantForGrowth(I,J,NZ,TFN6_vr,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
    RootPrimeAxsNum,TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4LeafPetolExpansion,TurgEff4CanopyResp)
  integer, intent(in) :: I,J,NZ
  REAL(R8), INTENT(OUT) :: TFN6_vr(JZ1)
  REAL(R8), INTENT(OUT) :: CNLFW,CPLFW
  real(r8), intent(out) :: CNSHW,CPSHW
  real(r8), intent(out) :: CNRTW,CPRTW
  real(r8), intent(out) :: RootPrimeAxsNum !Numer of primary root axes
  real(r8), intent(out) :: TFN5
  real(r8), intent(out) :: WaterStress4Groth
  real(r8), intent(out) :: Stomata_Stress
  real(r8), intent(out) :: TurgEff4LeafPetolExpansion,TurgEff4CanopyResp
  integer :: L,NR,N,NE
  real(r8) :: ACTVM,RTK,STK,TKCM,TKSM
  real(r8), parameter :: dscal=0.999992087_r8
  character(len=*), parameter :: subname='StagePlantForGrowth'
!     begin_execution

  associate(                                                               &
    rNCLeaf_pft                 => plt_allom%rNCLeaf_pft                  ,& !input  :maximum leaf N:C ratio, [g g-1]
    rNCSheath_pft               => plt_allom%rNCSheath_pft                ,& !input  :sheath N:C ratio, [g g-1]
    rPCLeaf_pft                 => plt_allom%rPCLeaf_pft                  ,& !input  :maximum leaf P:C ratio, [g g-1]
    rPCSheath_pft               => plt_allom%rPCSheath_pft                ,& !input  :sheath P:C ratio, [g g-1]
    CanopySapwoodC_pft          => plt_biom%CanopySapwoodC_pft            ,& !input  :canopy active stalk C, [g d-2]
    Myco_pft                    => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]
    MaxNumRootLays              => plt_site%MaxNumRootLays                ,& !input  :maximum root layer number,[-]
    NK                          => plt_site%NK                            ,& !input  :current hydrologically active layer, [-]
    NU                          => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    PSICanopyTurg_pft           => plt_ew%PSICanopyTurg_pft               ,& !input  :plant canopy turgor water potential, [MPa]
    PSICanopy_pft               => plt_ew%PSICanopy_pft                   ,& !input  :canopy total water potential, [Mpa]
    PlantPopulation_pft         => plt_site%PlantPopulation_pft           ,& !input  :plant population, [d-2]
    RCS_pft                     => plt_photo%RCS_pft                      ,& !input  :e-folding turgor pressure for stomatal resistance, [MPa]
    RootElms_pft                => plt_biom%RootElms_pft                  ,& !input  :plant root element mass, [g d-2]
    rNCRoot_pft                 => plt_allom%rNCRoot_pft                  ,& !input  :root N:C ratio, [gN gC-1]
    rPCRootr_pft                => plt_allom%rPCRootr_pft                 ,& !input  :root P:C ratio, [gP gC-1]
    StalkStrutElms_pft          => plt_biom%StalkStrutElms_pft            ,& !input  :canopy stalk structural element mass, [g d-2]
    TKC_pft                     => plt_ew%TKC_pft                         ,& !input  :canopy temperature, [K]
    TKS_vr                      => plt_ew%TKS_vr                          ,& !input  :mean annual soil temperature, [K]
    TempOffset_pft              => plt_pheno%TempOffset_pft               ,& !input  :adjustment of Arhhenius curves for plant thermal acclimation, [oC]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft        ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft    ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    rNCStalk_pft                => plt_allom%rNCStalk_pft                 ,& !input  :stalk N:C ratio, [gN gC-1]
    rPCStalk_pft                => plt_allom%rPCStalk_pft                 ,& !input  :stalk P:C ratio, [g g-1]
    k_fine_litr                 => pltpar%k_fine_litr                     ,& !input  :fine litter complex id
    k_woody_litr                => pltpar%k_woody_litr                    ,& !input  :woody litter complex id
    FracRootElmAlloc2Litr       => plt_allom%FracRootElmAlloc2Litr        ,& !inoput :C woody fraction in root,[-]
    FracWoodStalkElmAlloc2Litr  => plt_allom%FracWoodStalkElmAlloc2Litr   ,& !inoput :woody element allocation,[-]
    FracShootLeafAlloc2Litr     => plt_allom%FracShootLeafAlloc2Litr      ,& !inoput :woody element allocation, [-]
    FracShootPetolAlloc2Litr    => plt_allom%FracShootPetolAlloc2Litr     ,& !inoput :leaf element allocation,[-]
    RootBiomCPerPlant_pft       => plt_biom%RootBiomCPerPlant_pft         ,& !inoput :root C biomass per plant, [g p-1]
    CanopyLeafAreaZ_pft         => plt_morph%CanopyLeafAreaZ_pft          ,& !output :canopy layer leaf area, [m2 d-2]
    CanopyLeafCLyr_pft          => plt_biom%CanopyLeafCLyr_pft            ,& !output :canopy layer leaf C, [g d-2]
    CanopyStemAreaZ_pft         => plt_morph%CanopyStemAreaZ_pft          ,& !output :plant canopy layer stem area, [m2 d-2]
    Root1stXNumL_rpvr           => plt_morph%Root1stXNumL_rpvr            ,& !output :root layer number primary axes, [d-2]
    Root2ndXNumL_rpvr           => plt_morph%Root2ndXNumL_rpvr            ,& !output :root layer number axes, [d-2]
    RootCO2Autor_pvr            => plt_rbgc%RootCO2Autor_pvr              ,& !output :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr          => plt_rbgc%RootCO2EmisPot_pvr            ,& !output :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootN2Fix_pvr               => plt_bgcr%RootN2Fix_pvr                 ,& !output :root N2 fixation, [gN d-2 h-1]
    RootProteinC_pvr            => plt_biom%RootProteinC_pvr              ,& !output :root layer protein C, [gC d-2]
    RootRespPotent_pvr          => plt_rbgc%RootRespPotent_pvr             & !output :root respiration unconstrained by O2, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)
  D2: DO L=1,NumCanopyLayers1
    CanopyLeafAreaZ_pft(L,NZ) = 0._r8
    CanopyLeafCLyr_pft(L,NZ)  = 0._r8
    CanopyStemAreaZ_pft(L,NZ) = 0._r8
  ENDDO D2


  D6: DO L=1,NK
    D9: DO N=1,Myco_pft(NZ)    
      RootProteinC_pvr(N,L,NZ)   = 0._r8
      Root1stXNumL_rpvr(N,L,NZ)  = 0._r8
      Root2ndXNumL_rpvr(N,L,NZ)  = 0._r8
      RootRespPotent_pvr(N,L,NZ) = 0._r8
      RootCO2EmisPot_pvr(N,L,NZ) = 0._r8
      RootCO2Autor_pvr(N,L,NZ)   = 0._r8
    ENDDO D9
    RootN2Fix_pvr(L,NZ)          = 0._r8
  ENDDO D6

!
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     WTSTK,WVSTK=stalk,sapwood mass
!     FWOOD,FWODB=C woody fraction in stalk,other organs:0=woody,1=non-woody
!     CN*,CP*=N:C,P:C ratios in plant organs from PFT files
!     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!     *LF=leaf,*SHE=petiole,*STK=stalk,*RT=root
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!     FWOODN,FWOODP=N,P woody fraction in stalk:0=woody,1=non-woody
!
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))&
    .OR. StalkStrutElms_pft(ielmc,NZ).LE.ZERO4Groth_pft(NZ))THEN
    !not tree
    FracShootLeafAlloc2Litr(ielmc,k_fine_litr) = 1.0_r8
    FracWoodStalkElmAlloc2Litr(ielmc,k_fine_litr) = 1.0_r8
    FracRootElmAlloc2Litr(ielmc,k_fine_litr)      = 1.0_r8    
  ELSE
    !tree
    FracShootLeafAlloc2Litr(ielmc,k_fine_litr) = 1.0_r8
    FracWoodStalkElmAlloc2Litr(ielmc,k_fine_litr) = AMIN1(SQRT(CanopySapwoodC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ)),1._r8)
    FracRootElmAlloc2Litr(ielmc,k_fine_litr)      = AMIN1(SQRT(FRTX*CanopySapwoodC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ)),1._R8)
  ENDIF

  FracShootLeafAlloc2Litr(ielmc,k_woody_litr) = AZMAX1(1.0_r8-FracShootLeafAlloc2Litr(ielmc,k_fine_litr))
  FracWoodStalkElmAlloc2Litr(ielmc,k_woody_litr) = AZMAX1(1.0_r8-FracWoodStalkElmAlloc2Litr(ielmc,k_fine_litr))
  FracRootElmAlloc2Litr(ielmc,k_woody_litr)      = AZMAX1(1.0_r8-FracRootElmAlloc2Litr(ielmc,k_fine_litr))

  CNLFW=FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FracShootLeafAlloc2Litr(ielmc,k_fine_litr)*rNCLeaf_pft(NZ)
  CPLFW=FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FracShootLeafAlloc2Litr(ielmc,k_fine_litr)*rPCLeaf_pft(NZ)
  CNSHW=FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FracShootLeafAlloc2Litr(ielmc,k_fine_litr)*rNCSheath_pft(NZ)
  CPSHW=FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FracShootLeafAlloc2Litr(ielmc,k_fine_litr)*rPCSheath_pft(NZ)
  CNRTW=FracRootElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FracRootElmAlloc2Litr(ielmc,k_fine_litr)*rNCRoot_pft(NZ)
  CPRTW=FracRootElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FracRootElmAlloc2Litr(ielmc,k_fine_litr)*rPCRootr_pft(NZ)

  FracShootPetolAlloc2Litr(ielmc,1:NumOfPlantLitrCmplxs)=FracShootLeafAlloc2Litr(ielmc,1:NumOfPlantLitrCmplxs)

  FracShootLeafAlloc2Litr(ielmn,k_woody_litr) = AMIN1(FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNLFW,1._r8)
  FracShootLeafAlloc2Litr(ielmp,k_woody_litr) = AMIN1(FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPLFW,1._r8)
  
  FracShootPetolAlloc2Litr(ielmn,k_woody_litr)  = AMIN1(FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNSHW,1._r8)
  FracShootPetolAlloc2Litr(ielmp,k_woody_litr)  = AMIN1(FracShootLeafAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPSHW,1._r8)
  
  FracWoodStalkElmAlloc2Litr(ielmn,k_woody_litr)  = AMIN1(FracWoodStalkElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNRTW,1._r8)
  FracWoodStalkElmAlloc2Litr(ielmp,k_woody_litr)  = AMIN1(FracWoodStalkElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPRTW,1._r8)

  FracRootElmAlloc2Litr(ielmn,k_woody_litr)       = AMIN1(FracRootElmAlloc2Litr(ielmc,k_woody_litr)*rNCRoot_pft(NZ)/CNRTW,1._r8)
  FracRootElmAlloc2Litr(ielmp,k_woody_litr)       = AMIN1(FracRootElmAlloc2Litr(ielmc,k_woody_litr)*rPCRootr_pft(NZ)/CPRTW,1._r8)

  DO NE=2,NumPlantChemElms
    FracShootPetolAlloc2Litr(NE,k_fine_litr) = AZMAX1(1.0_r8-FracShootPetolAlloc2Litr(NE,k_woody_litr))
    FracShootLeafAlloc2Litr(NE,k_fine_litr)  = AZMAX1(1.0_r8-FracShootLeafAlloc2Litr(NE,k_woody_litr))
    FracWoodStalkElmAlloc2Litr(NE,k_fine_litr)  = AZMAX1(1.0_r8-FracWoodStalkElmAlloc2Litr(NE,k_woody_litr))
    FracRootElmAlloc2Litr(NE,k_fine_litr)       = AZMAX1(1.0_r8-FracRootElmAlloc2Litr(NE,k_woody_litr))
  ENDDO
!
!     SHOOT AND ROOT TEMPERATURE FUNCTIONS FOR MAINTENANCE
!     RESPIRATION FROM TEMPERATURES WITH OFFSETS FOR THERMAL ADAPTATION
!
!     TKC,TKCM=canopy temperature,canopy temp used in Arrhenius eqn
!     TKS_vr,TKSM=soil temperature,soil temp used in Arrhenius eqn
!     TempOffset_pft=shift in Arrhenius curve for thermal adaptation
!     TFN5,TFN6_vr=temperature function for canopy,root mntc respn (25 oC =1)
!     8.3143,710.0=gas constant,enthalpy
!     62500,195000,232500=energy of activn,high,low temp inactivn(KJ mol-1)
!
  TKCM=real_truncate(TKC_pft(NZ)+TempOffset_pft(NZ),1.e-3_r8)

  TFN5=calc_plant_maint_tempf(TKCM)  
  D7: DO L=NU,MaxNumRootLays
    TKSM       = real_truncate(TKS_vr(L)+TempOffset_pft(NZ),1.e-3_r8)
    TFN6_vr(L) = calc_plant_maint_tempf(TKSM)
  ENDDO D7
!
!     PRIMARY ROOT NUMBER
!
!     RootBiomCPerPlant_pft=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     RootPrimeAxsNum=multiplier for number of primary root axes
!
  RootBiomCPerPlant_pft(NZ) = AMAX1(dscal*RootBiomCPerPlant_pft(NZ),RootElms_pft(ielmc,NZ)/PlantPopulation_pft(NZ))
  RootPrimeAxsNum           = AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**0.667_r8)*PlantPopulation_pft(NZ)

  !
  !     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
  !     FROM CANOPY TURGOR
  !
  !     TurgEff4LeafPetolExpansion=turgor expansion,extension function
  !     PSICanopyTurg_pft,TurgPSIMin4OrganExtens=current,minimum canopy turgor potl for expansion,extension
  !     Stomata_Stress=stomatal resistance function of canopy turgor
  !     PSICanopy_pft=canopy water potential
  !     WaterStress4Groth=growth function of canopy water potential
  !     TurgEff4CanopyResp=expansion,extension function of canopy water potential
  !
  TurgEff4LeafPetolExpansion=real_truncate(AMIN1(1.0_r8,AZMAX1(PSICanopyTurg_pft(NZ)-TurgPSIMin4OrganExtens)),1.e-3_r8)

  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    !bryophyte, no turgor
    Stomata_Stress    = 1.0_r8
    WaterStress4Groth = EXP(0.05_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8))
  ELSE
    !others
    Stomata_Stress    = EXP(-PSICanopyTurg_pft(NZ)/RCS_pft(NZ))
    WaterStress4Groth = EXP(0.10_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8))
  ENDIF
  WaterStress4Groth               = real_truncate(WaterStress4Groth,1.e-3_r8)
  TurgEff4CanopyResp              = fRespWatSens(TurgEff4LeafPetolExpansion,iPlantRootProfile_pft(NZ))
  plt_biom%StomatalStress_pft(NZ) = Stomata_Stress

  call PrintInfo('end '//subname)
  end associate
  end subroutine StagePlantForGrowth

!----------------------------------------------------------------------------------------------------
  subroutine ComputeTotalBiom(I,J,NZ)

  integer, intent(in) :: I,J,NZ
  
  integer :: L,N
!     begin_execution
  associate(                                                    &
    Myco_pft               => plt_morph%Myco_pft               ,& !input  :mycorrhizal type (no or yes),[-]
    MaxSoiL4Root_pft       => plt_morph%MaxSoiL4Root_pft       ,& !input  :maximum soil layer number for all root axes,[-]
    NU                     => plt_site%NU                      ,& !input  :current soil surface layer number, [-]
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr  ,& !input  :root layer nonstructural element, [g d-2]
    PopuRootMycoC_pvr      => plt_biom% PopuRootMycoC_pvr       & !inoput :root layer C, [gC d-2]
  )
  !prepare for disturbance
  CALL SumPlantBiomStates(I,J,NZ,'computotb')

  D3451: DO N=1,Myco_pft(NZ)
    DO  L=NU,MaxSoiL4Root_pft(NZ)
      PopuRootMycoC_pvr(N,L,NZ)=PopuRootMycoC_pvr(N,L,NZ)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
    enddo
  ENDDO D3451

  end associate

  end subroutine ComputeTotalBiom

!----------------------------------------------------------------------------------------------------
  subroutine AccumulateStates(I,J,NZ)
  
  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L,NR,N,NE,NB

!     begin_execution
  associate(                                                          &
    ZEROS                     => plt_site%ZEROS                      ,& !input  :threshold zero for numerical stability,[-]  
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch  ,& !input  :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch  ,& !input  :branch nodule structural element, [g d-2]
    CanopyStalkArea_lbrch     => plt_morph%CanopyStalkArea_lbrch     ,& !input  :plant canopy layer branch stem area, [m2 d-2]
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch     ,& !input  :plant branch leaf + sheath C, [g d-2]
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft          ,& !input  :maximum soil layer number for all root axes,[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    RootH2PO4Uptake_pft       => plt_rbgc%RootH2PO4Uptake_pft        ,& !input  :total root uptake of PO4, [g d-2 h-1]
    RootHPO4Uptake_pft        => plt_rbgc%RootHPO4Uptake_pft         ,& !input  :total root uptake of HPO4, [g d-2 h-1]
    RootMycoExudElms_pft      => plt_rbgc%RootMycoExudElms_pft       ,& !input  :total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
    RootN2Fix_pft             => plt_rbgc%RootN2Fix_pft              ,& !input  :total root N2 fixation, [g d-2 h-1]
    RootNH4Uptake_pft         => plt_rbgc%RootNH4Uptake_pft          ,& !input  :total root uptake of NH4, [g d-2 h-1]
    RootNO3Uptake_pft         => plt_rbgc%RootNO3Uptake_pft          ,& !input  :total root uptake of NO3, [g d-2 h-1]
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr    ,& !input  :root layer nonstructural element, [g d-2]
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr    ,& !input  :root layer nodule element, [g d-2]
    CanopyN2Fix_pft           => plt_rbgc%CanopyN2Fix_pft            ,& !input :total canopy N2 fixation, [g d-2 h-1]                    
    SeedSitesSet_brch         => plt_morph%SeedSitesSet_brch         ,& !input  :branch grain number, [d-2]
    SapwoodBiomassC_brch      => plt_biom%SapwoodBiomassC_brch       ,& !input  :branch live stalk C, [gC d-2]
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    CanopyNodulNonstElms_pft  => plt_biom%CanopyNodulNonstElms_pft   ,& !inoput :canopy nodule nonstructural element, [g d-2]
    CanopyStemAreaZ_pft       => plt_morph%CanopyStemAreaZ_pft       ,& !inoput :plant canopy layer stem area, [m2 d-2]
    PlantExudElm_CumYr_pft    => plt_rbgc%PlantExudElm_CumYr_pft     ,& !inoput :total net root element uptake (+ve) - exudation (-ve), [gC d-2 ]
    PlantN2Fix_CumYr_pft      => plt_bgcr%PlantN2Fix_CumYr_pft       ,& !inoput :total plant N2 fixation, [g d-2 ]
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft   ,& !inoput :net root element uptake (+ve) - exudation (-ve), [gC d-2 h-1]
    RootUptk_N_CumYr_pft      => plt_rbgc%RootUptk_N_CumYr_pft       ,& !inoput :cumulative plant N uptake, [gN d-2]
    RootUptk_P_CumYr_pft      => plt_rbgc%RootUptk_P_CumYr_pft       ,& !inoput :cumulative plant P uptake, [gP d-2]
    VcMaxRubiscoRef_brch      => plt_photo%VcMaxRubiscoRef_brch      ,& !input  :maximum rubisco carboxylation rate at reference temperature, [umol g-1 h-1]    
    VoMaxRubiscoRef_brch      => plt_photo%VoMaxRubiscoRef_brch      ,& !input  :maximum rubisco oxygenation rate at reference temperature, [umol g-1 h-1]    
    VcMaxPEPCarboxyRef_brch   => plt_photo%VcMaxPEPCarboxyRef_brch   ,& !input  :reference maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]        
    fNCLFW_brch               => plt_pheno%fNCLFW_brch               ,& !input : NC ratio of growing leaf on branch, [gN/gC]
    fPCLFW_brch               => plt_pheno%fPCLFW_brch               ,& !input : PC ratio of growing leaf on branch, [gP/gC]
    LeafAreaLive_brch         => plt_morph%LeafAreaLive_brch         ,& !output :branch leaf area, [m2 d-2]
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch         ,& !output :branch leaf structural element mass, [g d-2]
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch       ,& !output :branch sheath structural element, [g d-2]
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft        ,& !output :plant canopy leaf area, [m2 d-2]
    CanopyLeafSheathC_pft     => plt_biom%CanopyLeafSheathC_pft      ,& !output :canopy leaf + sheath C, [g d-2]
    CanopySeedNum_pft         => plt_morph%CanopySeedNum_pft         ,& !output :canopy grain number, [d-2]
    CanopySapwoodC_pft        => plt_biom%CanopySapwoodC_pft         ,& !output :canopy active stalk C, [g d-2]
    CanopyStemArea_pft        => plt_morph%CanopyStemArea_pft        ,& !output :plant stem area, [m2 d-2]
    EarStrutElms_pft          => plt_biom%EarStrutElms_pft           ,& !output :canopy ear structural element, [g d-2]
    GrainStrutElms_pft        => plt_biom%GrainStrutElms_pft         ,& !output :canopy grain structural element, [g d-2]
    HuskStrutElms_pft         => plt_biom%HuskStrutElms_pft          ,& !output :canopy husk structural element mass, [g d-2]
    LeafStrutElms_pft         => plt_biom%LeafStrutElms_pft          ,& !output :canopy leaf structural element mass, [g d-2]
    PetoleStrutElms_pft       => plt_biom%PetoleStrutElms_pft        ,& !output :canopy sheath structural element mass, [g d-2]
    StalkRsrvElms_pft         => plt_biom%StalkRsrvElms_pft          ,& !output :canopy reserve element mass, [g d-2]
    StalkStrutElms_pft        => plt_biom%StalkStrutElms_pft         ,& !output :canopy stalk structural element mass, [g d-2]
    fNCLFW_pft                => plt_pheno%fNCLFW_pft                ,& !output : NC ratio of growing leaf, [gN/gC]
    fPCLFW_pft                => plt_pheno%fPCLFW_pft                 & !output : PC ratio of growing leaf, [gP/gC]

  )
!
!     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
!
!
  DO NB=1,NumOfBranches_pft(NZ)    
    if(isclose(CanopyLeafSheathC_brch(NB,NZ),0._r8))then
      LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
      PetoleStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
      LeafAreaLive_brch(NB,NZ)                     = 0._r8
    endif
 ENDDO  

  CanopySapwoodC_pft(NZ)   = 0._r8
  CanopyLeafSheathC_pft(NZ) = 0._r8
  CanopySeedNum_pft(NZ)    = 0._r8
  CanopyLeafArea_pft(NZ)   = 0._r8
  CanopyStemArea_pft(NZ)   = 0._r8

  DO NB=1,NumOfBranches_pft(NZ)        
    CanopySapwoodC_pft(NZ)     = CanopySapwoodC_pft(NZ)+SapwoodBiomassC_brch(NB,NZ)
    CanopyLeafSheathC_pft(NZ)   = CanopyLeafSheathC_pft(NZ) +CanopyLeafSheathC_brch(NB,NZ)
    CanopySeedNum_pft(NZ)      = CanopySeedNum_pft(NZ)+SeedSitesSet_brch(NB,NZ)
    CanopyLeafArea_pft(NZ)     = CanopyLeafArea_pft(NZ)+LeafAreaLive_brch(NB,NZ)
    fNCLFW_pft(NZ)             = fNCLFW_pft(NZ)+fNCLFW_brch(NB,NZ)*LeafAreaLive_brch(NB,NZ)
    fPCLFW_pft(NZ)             = fPCLFW_pft(NZ)+fPCLFW_brch(NB,NZ)*LeafAreaLive_brch(NB,NZ)
    DO L=1,NumCanopyLayers1
      CanopyStemAreaZ_pft(L,NZ)=CanopyStemAreaZ_pft(L,NZ)+CanopyStalkArea_lbrch(L,NB,NZ)
    ENDDO
  ENDDO
  if(CanopyLeafArea_pft(NZ).GT.ZEROs)then
    fNCLFW_pft(NZ)=fNCLFW_pft(NZ)/CanopyLeafArea_pft(NZ)
    fPCLFW_pft(NZ)=fPCLFW_pft(NZ)/CanopyLeafArea_pft(NZ)
  endif

!
!     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
!
!     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
!
!     ACCUMULATE TOTAL SOIL-PLANT C,N,P EXCHANGE
!
!     PlantN2Fix_CumYr_pft=cumulative PFT N2 fixation
!     PlantRootSoilElmNetX_pft= > 0 taking from soil 
  PlantRootSoilElmNetX_pft(1:NumPlantChemElms,NZ)=RootMycoExudElms_pft(1:NumPlantChemElms,NZ)

  PlantRootSoilElmNetX_pft(ielmn,NZ)=PlantRootSoilElmNetX_pft(ielmn,NZ)+&
    RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ)+RootN2Fix_pft(NZ)
  PlantRootSoilElmNetX_pft(ielmp,NZ)=PlantRootSoilElmNetX_pft(ielmp,NZ)+&
    RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ)

  PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)=PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)+&
    RootMycoExudElms_pft(1:NumPlantChemElms,NZ)

  RootUptk_N_CumYr_pft(NZ)=RootUptk_N_CumYr_pft(NZ)+RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ) + &
    RootMycoExudElms_pft(ielmn,NZ)
  RootUptk_P_CumYr_pft(NZ)=RootUptk_P_CumYr_pft(NZ)+RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ) + &
    RootMycoExudElms_pft(ielmp,NZ)

  PlantN2Fix_CumYr_pft(NZ)=PlantN2Fix_CumYr_pft(NZ)+RootN2Fix_pft(NZ)+CanopyN2Fix_pft(NZ)
  end associate
  end subroutine AccumulateStates

!----------------------------------------------------------------------------------------------------
  subroutine ZeroPlantStates(I,J,NZ)  
  !
  !for annual plants
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NE
  DO NE=1,NumPlantChemElms
    plt_biom%CanopyNonstElms_pft(NE,NZ)=0._r8
    plt_biom%StalkRsrvElms_pft(NE,NZ)=0._r8
    plt_biom%ShootElms_pft(NE,NZ)=0._r8
    plt_biom%PetoleStrutElms_pft(NE,NZ)=0._r8
    plt_biom%StalkStrutElms_pft(NE,NZ)=0._r8
    plt_biom%LeafStrutElms_pft(NE,NZ)=0._r8
    plt_biom%HuskStrutElms_pft(NE,NZ)=0._r8
    plt_biom%GrainStrutElms_pft(NE,NZ)=0._r8
    plt_biom%EarStrutElms_pft(NE,NZ)=0._r8
    plt_biom%ShootNoduleElms_pft(NE,NZ)=0._r8
    plt_biom%RootNoduleElms_pft(NE,NZ)=0._r8
    plt_rbgc%PlantRootSoilElmNetX_pft(NE,NZ)=0._r8
  ENDDO

  plt_biom%CanopySapwoodC_pft(NZ)=0._r8
  plt_biom%CanopyLeafSheathC_pft(NZ) =0._r8
  plt_morph%CanopySeedNum_pft(NZ) =0._r8
  plt_morph%CanopyLeafArea_pft(NZ)=0._r8
  plt_morph%CanopyStemArea_pft(NZ)=0._r8
  plt_morph%CanopyStemAreaZ_pft(1:NumCanopyLayers1,NZ)=0._r8

  end subroutine ZeroPlantStates

  ![tail]
end module grosubsMod
