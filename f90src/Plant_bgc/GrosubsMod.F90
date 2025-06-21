module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod, only : safe_adb,AZMAX1
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod,       only: lverb
  use EcoSiMParDataMod,    only: pltpar
  use RootMod,             only: RootBGCModel
  use PlantNonstElmDynMod, only: PlantNonstElmTransfer
  use PlantDebugMod,       only: PrintRootTracer
  use EcosimConst
  use GrosubPars
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
  public :: InitGrosub
  contains

  subroutine InitGrosub(NumGrowthStages,MaxNumRootAxes)

  implicit none
  integer, intent(out) :: NumGrowthStages,MaxNumRootAxes

  call InitVegPars(pltpar)
  
  NumGrowthStages = pltpar%NumGrowthStages
  MaxNumRootAxes  = pltpar%MaxNumRootAxes

  end subroutine InitGrosub
!------------------------------------------------------------------------------------------

  subroutine GrowPlants(I,J)
!
!     THIS subroutine CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
!
  use PlantDisturbsMod, only : RemoveBiomassByDisturbance
  implicit none
  integer, intent(in) :: I, J

  real(r8) :: CanopyHeight_copy(JP1)
  integer :: L,K,M
  integer :: NZ,NE
! begin_execution
  associate(                                          &
    IsPlantActive_pft => plt_pheno%IsPlantActive_pft, &  !input
    NP                => plt_site%NP,                 &  !input
    NP0               => plt_site%NP0,                &  !input
    CanopyHeight_pft  => plt_morph%CanopyHeight_pft   &  !inoput
  )
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
!
!     INITIALIZE SENESCENCE ARRAYS
  DO NZ=1,NP0
    CanopyHeight_copy(NZ)            = CanopyHeight_pft(NZ)
    CanopyHeight_pft(NZ)             = 0._r8
    plt_biom%RootMassElm_pvr(:,:,NZ) = 0._r8
    plt_rbgc%canopy_growth_pft(NZ)   = 0._r8
  ENDDO  
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
  D9985: DO NZ=1,NP

    ! IsPlantActive_pft= flag for living pft
    IF(IsPlantActive_pft(NZ).EQ.iActive)THEN
      if(lverb)write(*,*)'GrowPlant'
      call GrowOnePlant(I,J,NZ,CanopyHeight_copy)
    else
      call AccumulateStates(I,J,NZ)
    ENDIF

    if(lverb)write(*,*)'HARVEST STANDING DEAD'
    call RemoveBiomassByDisturbance(I,J,NZ)

    call PrintRootTracer(I,J,NZ,'afdist')
  ENDDO D9985
!
  if(lverb)write(*,*)'TRANSFORMATIONS IN LIVING OR DEAD PLANT POPULATIONS'
  call LiveDeadTransformation(I,J)
  
  DO NZ=1,NP
    call SumPlantBiom(I,J,NZ,'exgrosubs')
    call PrintRootTracer(I,J,NZ,'afddeadt')    
  ENDDO
  end associate
  END subroutine GrowPlants

!------------------------------------------------------------------------------------------

  subroutine LiveDeadTransformation(I,J)
  implicit none
  integer, intent(in) :: I,J

  integer :: L,K,NZ,M,NE,NB
  real(r8) :: XFRC,XFRN,XFRP,XFRE
!     begin_execution

  associate(                                                                   &
    CH4ByFire_CumYr_pft            => plt_distb%CH4ByFire_CumYr_pft,           &  !input
    CO2ByFire_CumYr_pft            => plt_distb%CO2ByFire_CumYr_pft,           &  !input
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr,            &  !input
    EcoHavstElmntCum_pft           => plt_distb%EcoHavstElmntCum_pft,          &  !input
    GrossCO2Fix_pft                => plt_bgcr%GrossCO2Fix_pft,                &  !input
    GrossResp_pft                  => plt_bgcr%GrossResp_pft,                  &  !input
    HourReq4LeafOut_brch           => plt_pheno%HourReq4LeafOut_brch,          &  !input
    Hours4Leafout_brch             => plt_pheno%Hours4Leafout_brch,            &  !input
    IsPlantActive_pft              => plt_pheno%IsPlantActive_pft,             &  !input
    MaxNumRootLays                 => plt_site%MaxNumRootLays,                 &  !input
    N2ObyFire_CumYr_pft            => plt_distb%N2ObyFire_CumYr_pft,           &  !input
    NH3Emis_CumYr_pft              => plt_bgcr%NH3Emis_CumYr_pft,              &  !input
    NH3byFire_CumYr_pft            => plt_distb%NH3byFire_CumYr_pft,           &  !input
    NP0                            => plt_site%NP0,                            &  !input
    NetCumElmntFlx2Plant_pft       => plt_pheno%NetCumElmntFlx2Plant_pft,      &  !input
    NodulStrutElms_pft             => plt_biom%NodulStrutElms_pft,             &  !input
    NumOfBranches_pft              => plt_morph%NumOfBranches_pft,             &  !input
    PO4byFire_CumYr_pft            => plt_distb%PO4byFire_CumYr_pft,           &  !input
    PlantExudElm_CumYr_pft         => plt_rbgc%PlantExudElm_CumYr_pft,         &  !input
    PlantN2Fix_CumYr_pft           => plt_bgcr%PlantN2Fix_CumYr_pft,           &  !input
    RootElms_pft                   => plt_biom%RootElms_pft,                   &  !input
    SeasonalNonstElms_pft          => plt_biom%SeasonalNonstElms_pft,          &  !input
    ShootStrutElms_pft             => plt_biom%ShootStrutElms_pft,             &  !input
    StandDeadStrutElms_pft         => plt_biom%StandDeadStrutElms_pft,         &  !input
    doPlantLeafOut_brch            => plt_pheno%doPlantLeafOut_brch,           &  !input
    fTCanopyGroth_pft              => plt_pheno%fTCanopyGroth_pft,             &  !input
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft,         &  !input
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft,     &  !input
    iYearCurrent                   => plt_site%iYearCurrent,                   &  !input
    ElmBalanceCum_pft              => plt_site%ElmBalanceCum_pft,              &  !inoput
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft,     &  !inoput
    LitrfalStrutElms_pft           => plt_bgcr%LitrfalStrutElms_pft,           &  !inoput
    LitrfalStrutElms_pvr           => plt_bgcr%LitrfalStrutElms_pvr,           &  !inoput
    NetPrimProduct_pft             => plt_bgcr%NetPrimProduct_pft,             &  !inoput
    StandDeadKCompElms_pft         => plt_biom%StandDeadKCompElms_pft,         &  !inoput
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft, &  !inoput
    doInitPlant_pft                => plt_pheno%doInitPlant_pft,               &  !inoput
    k_fine_litr                    => pltpar%k_fine_litr,                      &  !inoput
    k_woody_litr                   => pltpar%k_woody_litr,                     &  !inoput
    PlantinDepz_pft                => plt_morph%PlantinDepz_pft,               &  !output
    iDayPlanting_pft               => plt_distb%iDayPlanting_pft,              &  !output
    iYearPlanting_pft              => plt_distb%iYearPlanting_pft              &  !output
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
!     LitrFall FROM STANDING DEAD
!
!     XFRC,XFRN,XFRP=LitrFall from standing dead
!     fTCanopyGroth_pft=temperature function for canopy growth
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     CSNC,ZSNC,PSNC=C,N,P LitrFall
!
    
    D6235: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        XFRE=1.5814E-05_r8*fTCanopyGroth_pft(NZ)*StandDeadKCompElms_pft(NE,M,NZ)
        IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.iPlantRootProfile_pft(NZ).LE.1)THEN
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XFRE
        ELSE
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XFRE
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
          SurfLitrfalStrutElms_CumYr_pft(NE,NZ)=SurfLitrfalStrutElms_CumYr_pft(NE,NZ)+LitrfalStrutElms_pvr(NE,M,K,0,NZ)
        ENDDO
      ENDDO D6431  
    ENDDO
!   the following purposely starts from layer 0.
    LitrfalStrutElms_pft(1:NumPlantChemElms,NZ)=0._r8
    D8955: DO L=0,MaxNumRootLays
      DO K=1,pltpar%NumOfPlantLitrCmplxs      
        D6430: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pft(NE,NZ)=LitrfalStrutElms_pft(NE,NZ)+LitrfalStrutElms_pvr(NE,M,K,L,NZ)
          enddo         
        ENDDO D6430      
      ENDDO
    ENDDO D8955    

    DO NE=1,NumPlantChemElms
      LitrfalStrutElms_CumYr_pft(NE,NZ)=LitrfalStrutElms_CumYr_pft(NE,NZ)+LitrfalStrutElms_pft(NE,NZ)
    ENDDO

!
!     PLANT C BALANCE = TOTAL C STATE VARIABLES + TOTAL
!     AUTOTROPHIC RESPIRATION + TOTAL LitrFall - TOTAL EXUDATION
!     - TOTAL CO2 FIXATION
!
!     BALC=PFT C balance, dplant/dt=input-output
!     starting from time zero, plant(t)-plant(t0)=input-output
!     input: net carbon fixation, 
!     output: litterfal, exudation, harvest, mortality
!
!     WTSHT,WTRT,WTND,WTRVC,WTSTG=PFT shoot,root,bacteria,storage,standing dead C
!     NetPrimProduct_pft=cumulative PFT NPP
!     TCSNC=cumulative PFT C LitrFall
!     TCUPTK=cumulative PFT root-soil C exchange
!     RSETC=cumulative C balance from previous year
!     THVSTC=cumulative PFT C removed from ecosystem from previous year
!     HVSTC=total PFT C removed from ecosystem in current year
!     CO2ByFire_CumYr_pft,CH4ByFire_CumYr_pft=CO2,CH4 emission from disturbance
!     LitrfalStrutElms_CumYr_pft= >0 to the environment (soil + surface)
!   GrossResp_pft < 0 respired into atmosphere
    NetPrimProduct_pft(NZ) = GrossCO2Fix_pft(NZ)+GrossResp_pft(NZ)

    IF(IsPlantActive_pft(NZ).EQ.iActive)THEN
      !check for living plant
      DO NE=1,NumPlantChemElms
        ElmBalanceCum_pft(NE,NZ)=ShootStrutElms_pft(NE,NZ)+RootElms_pft(NE,NZ)+NodulStrutElms_pft(NE,NZ) &
          +SeasonalNonstElms_pft(NE,NZ) +StandDeadStrutElms_pft(NE,NZ)     &   !add biomass by components  
          -LitrfalStrutElms_CumYr_pft(NE,NZ)-PlantExudElm_CumYr_pft(NE,NZ) &   !add fluxes
          -NetCumElmntFlx2Plant_pft(NE,NZ)+EcoHavstElmntCum_pft(NE,NZ)
      ENDDO
      !add more fluxes
      ElmBalanceCum_pft(ielmc,NZ) = ElmBalanceCum_pft(ielmc,NZ)-NetPrimProduct_pft(NZ) &
        -CO2ByFire_CumYr_pft(NZ)-CH4ByFire_CumYr_pft(NZ)
!
!     PLANT N BALANCE = TOTAL N STATE VARIABLES + TOTAL N LitrFall
!     - TOTAL N UPTAKE FROM SOIL - TOTAL N ABSORPTION FROM ATMOSPHERE
!
!     BALN=PFT N balance
!     WTSHN,WTRTN,WTNDN,WTRVN,WTSTGN=PFT shoot,root,bacteria,storage,standing dead N
!     TZSNC=cumulative PFT N LitrFall
!     TZUPTK=cumulative PFT root-soil N exchange
!     NH3Emis_CumYr_pft=cumulative NH3 exchange
!     RSETN=cumulative N balance from previous year
!     THVSTN=cumulative PFT N removed from ecosystem from previous year
!     HVSTN=total PFT N removed from ecosystem in current year
!     NH3byFire_CumYr_pft,N2ObyFire_CumYr_pft=NH3,N2O emission from disturbance
!     PlantN2Fix_CumYr_pft=cumulative PFT N2 fixation
!
      ElmBalanceCum_pft(ielmn,NZ)=ElmBalanceCum_pft(ielmn,NZ)-NH3Emis_CumYr_pft(NZ) &
        -NH3byFire_CumYr_pft(NZ)-N2ObyFire_CumYr_pft(NZ)-PlantN2Fix_CumYr_pft(NZ)
!
!     PLANT P BALANCE = TOTAL P STATE VARIABLES + TOTAL P LitrFall
!     - TOTAL P UPTAKE FROM SOIL
!
!     BALP=PFT N balance
!     WTSHP,WTRTP,WTNDP,WTRVP,WTSTGP=PFT shoot,root,bacteria,storage,standing dead P
!     TPSNC=cumulative PFT P LitrFall
!     TPUPTK=cumulative PFT root-soil P exchange
!     RSETP=cumulative P balance from previous year
!     THVSTP=cumulative PFT P removed from ecosystem from previous year
!     HVSTP=total PFT P removed from ecosystem in current year
!     PO4byFire_CumYr_pft=PO4 emission from disturbance
!
      ElmBalanceCum_pft(ielmp,NZ)=ElmBalanceCum_pft(ielmp,NZ)-PO4byFire_CumYr_pft(NZ)
    ENDIF
  ENDDO D9975
  end associate
  end subroutine LiveDeadTransformation
!------------------------------------------------------------------------------------------

  subroutine GrowOnePlant(I,J,NZ,CanopyHeight_copy)
  !
  !Description
  !plant growth
  use PlantDisturbsMod, only : RemoveBiomByMgmt
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)

  real(r8)  :: CanopyN2Fix_pft(JP1)
  integer  :: NB
  integer  :: BegRemoblize
  real(r8) :: TFN6_vr(JZ1)
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT
  real(r8) :: RootPrimeAxsNum
  real(r8) :: TFN5
  real(r8) :: WaterStress4Groth
  real(r8) :: Stomata_Stress
  real(r8) :: WFNS,CanTurgPSIFun4Expans
  real(r8) :: RootSinkC_vr(jroots,JZ1),RootSinkC(jroots)
  real(r8) :: Root1stSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: Root2ndSink_pvr(jroots,JZ1,pltpar%MaxNumRootAxes)
  
! begin_execution
  associate(                                               &
    NumOfBranches_pft    => plt_morph%NumOfBranches_pft,   &  !input
    iPlantRootState_pft  => plt_pheno%iPlantRootState_pft, &  !input
    iPlantShootState_pft => plt_pheno%iPlantShootState_pft &  !input
  )

  IF(iPlantShootState_pft(NZ).EQ.iLive .OR. iPlantRootState_pft(NZ).EQ.iLive)THEN
    CanopyN2Fix_pft(NZ) = 0._r8
    BegRemoblize        = 0
    
    call StagePlantForGrowth(I,J,NZ,TFN6_vr,CNLFW,CPLFW,&
      CNSHW,CPSHW,CNRTW,CPRTW,RootPrimeAxsNum,TFN5,WaterStress4Groth,Stomata_Stress,WFNS,CanTurgPSIFun4Expans)
!
!     CALCULATE GROWTH OF EACH BRANCH

    DO  NB=1,NumOfBranches_pft(NZ)
      if(lverb)write(*,*)'GrowOneBranch'
      call GrowOneBranch(I,J,NB,NZ,TFN6_vr,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WaterStress4Groth,Stomata_Stress,WFNS,CanTurgPSIFun4Expans,PTRT,CanopyN2Fix_pft,BegRemoblize)
    ENDDO

    call PrintRootTracer(I,J,NZ,'afgrowbranch')
    if(lverb)write(*,*)'RootBGCModel'
    call RootBGCModel(I,J,NZ,BegRemoblize,PTRT,TFN6_vr,CNRTW,CPRTW,RootPrimeAxsNum, &
      RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)

    call PrintRootTracer(I,J,NZ,'afroot')
    if(lverb)write(*,*)'PlantNonstElmTransfer'
    call PlantNonstElmTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC,BegRemoblize)

    if(lverb)write(*,*)'ComputeTotalBiom'
    call ComputeTotalBiom(I,J,NZ)
  ENDIF
!
  if(lverb)write(*,*)'RemoveBiomByMgmt'
  call RemoveBiomByMgmt(I,J,NZ)

  call PrintRootTracer(I,J,NZ,'afbiomrm')
!
!     RESET DEAD BRANCHES
  call ResetDeadPlant(I,J,NZ)

  call AccumulateStates(I,J,NZ,CanopyN2Fix_pft)

  call PrintRootTracer(I,J,NZ,'afaccum')
  end associate
  end subroutine GrowOnePlant

!------------------------------------------------------------------------------------------

  subroutine StagePlantForGrowth(I,J,NZ,TFN6_vr,CNLFW,CPLFW,CNSHW,&
    CPSHW,CNRTW,CPRTW,RootPrimeAxsNum,TFN5,WaterStress4Groth,Stomata_Stress,WFNS,CanTurgPSIFun4Expans)
  integer, intent(in) :: I,J,NZ
  REAL(R8), INTENT(OUT):: TFN6_vr(JZ1)
  REAL(R8), INTENT(OUT) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8), intent(out) :: RootPrimeAxsNum
  real(r8), intent(out) :: TFN5
  real(r8), intent(out) :: WaterStress4Groth
  real(r8), intent(out) :: Stomata_Stress
  real(r8), intent(out) :: WFNS,CanTurgPSIFun4Expans
  integer :: L,NR,N,NE
  real(r8) :: ACTVM,RTK,STK,TKCM,TKSM
!     begin_execution

  associate(                                                              &
    CNLF_pft                    => plt_allom%CNLF_pft,                    &  !input
    CNSHE_pft                   => plt_allom%CNSHE_pft,                   &  !input
    CPLF_pft                    => plt_allom%CPLF_pft,                    &  !input
    CPSHE_pft                   => plt_allom%CPSHE_pft,                   &  !input
    CanopyStalkC_pft            => plt_biom%CanopyStalkC_pft,             &  !input
    Myco_pft                      => plt_morph%Myco_pft,                      &  !input
    MaxNumRootLays              => plt_site%MaxNumRootLays,               &  !input
    NK                          => plt_site%NK,                           &  !input
    NU                          => plt_site%NU,                           &  !input
    PSICanopyTurg_pft           => plt_ew%PSICanopyTurg_pft,              &  !input
    PSICanopy_pft               => plt_ew%PSICanopy_pft,                  &  !input
    PlantPopulation_pft         => plt_site%PlantPopulation_pft,          &  !input
    RCS_pft                     => plt_photo%RCS_pft,                     &  !input
    RootElms_pft                => plt_biom%RootElms_pft,                 &  !input
    RootrNC_pft                 => plt_allom%RootrNC_pft,                 &  !input
    RootrPC_pft                 => plt_allom%RootrPC_pft,                 &  !input
    StalkStrutElms_pft          => plt_biom%StalkStrutElms_pft,           &  !input
    TKC_pft                     => plt_ew%TKC_pft,                        &  !input
    TKS_vr                      => plt_ew%TKS_vr,                         &  !input
    TempOffset_pft              => plt_pheno%TempOffset_pft,              &  !input
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft,               &  !input
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft,       &  !input
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft,   &  !input
    rNCStalk_pft                => plt_allom%rNCStalk_pft,                &  !input
    rPCStalk_pft                => plt_allom%rPCStalk_pft,                &  !input
    FracRootElmAlloc2Litr       => plt_allom%FracRootElmAlloc2Litr,       &  !inoput
    FracRootStalkElmAlloc2Litr  => plt_allom%FracRootStalkElmAlloc2Litr,  &  !inoput
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr,  &  !inoput
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr, &  !inoput
    RootBiomCPerPlant_pft       => plt_biom%RootBiomCPerPlant_pft,        &  !inoput
    k_fine_litr                 => pltpar%k_fine_litr,                    &  !inoput
    k_woody_litr                => pltpar%k_woody_litr,                   &  !inoput
    CanopyLeafAreaZ_pft         => plt_morph%CanopyLeafAreaZ_pft,         &  !output
    CanopyLeafCLyr_pft          => plt_biom%CanopyLeafCLyr_pft,           &  !output
    CanopyStemAreaZ_pft         => plt_morph%CanopyStemAreaZ_pft,         &  !output
    Root1stXNumL_pvr            => plt_morph%Root1stXNumL_pvr,            &  !output
    Root2ndXNum_pvr             => plt_morph%Root2ndXNum_pvr,             &  !output
    RootCO2Autor_pvr            => plt_rbgc%RootCO2Autor_pvr,             &  !output
    RootCO2EmisPot_pvr          => plt_rbgc%RootCO2EmisPot_pvr,           &  !output
    RootN2Fix_pvr               => plt_bgcr%RootN2Fix_pvr,                &  !output
    RootProteinC_pvr            => plt_biom%RootProteinC_pvr,             &  !output
    RootRespPotent_pvr          => plt_rbgc%RootRespPotent_pvr            &  !output
   )

  D2: DO L=1,NumOfCanopyLayers1
    CanopyLeafAreaZ_pft(L,NZ)=0._r8
    CanopyLeafCLyr_pft(L,NZ)=0._r8
    CanopyStemAreaZ_pft(L,NZ)=0._r8
  ENDDO D2


  D6: DO L=1,NK
    D9: DO N=1,Myco_pft(NZ)    
      RootProteinC_pvr(N,L,NZ)   = 0._r8
      Root1stXNumL_pvr(N,L,NZ)   = 0._r8
      Root2ndXNum_pvr(N,L,NZ)    = 0._r8
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
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))&
    .OR.StalkStrutElms_pft(ielmc,NZ).LE.ZERO4Groth_pft(NZ))THEN
    FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr) = 1.0_r8
    FracRootStalkElmAlloc2Litr(ielmc,k_fine_litr) = 1.0_r8
    FracRootElmAlloc2Litr(ielmc,k_fine_litr)      = 1.0_r8
  ELSE
    FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr) = 1.0_r8
    FracRootStalkElmAlloc2Litr(ielmc,k_fine_litr) = SQRT(CanopyStalkC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ))
    FracRootElmAlloc2Litr(ielmc,k_fine_litr)      = SQRT(FRTX*CanopyStalkC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ))
  ENDIF

  FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr) = 1.0_r8-FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr)
  FracRootStalkElmAlloc2Litr(ielmc,k_woody_litr) = 1.0_r8-FracRootStalkElmAlloc2Litr(ielmc,k_fine_litr)
  FracRootElmAlloc2Litr(ielmc,k_woody_litr)      = 1.0_r8-FracRootElmAlloc2Litr(ielmc,k_fine_litr)

  CNLFW=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr)*CNLF_pft(NZ)
  CPLFW=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr)*CPLF_pft(NZ)

  CNSHW=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr)*CNSHE_pft(NZ)
  CPSHW=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr)*CPSHE_pft(NZ)

  CNRTW=FracRootElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FracRootElmAlloc2Litr(ielmc,k_fine_litr)*RootrNC_pft(NZ)
  CPRTW=FracRootElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FracRootElmAlloc2Litr(ielmc,k_fine_litr)*RootrPC_pft(NZ)

  FracShootStalkElmAlloc2Litr(ielmc,1:NumOfPlantLitrCmplxs)=FracShootLeafElmAlloc2Litr(ielmc,1:NumOfPlantLitrCmplxs)

  FracShootStalkElmAlloc2Litr(ielmn,k_woody_litr)=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNLFW
  FracShootStalkElmAlloc2Litr(ielmp,k_woody_litr)=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPLFW

  FracShootLeafElmAlloc2Litr(ielmn,k_woody_litr)=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNSHW
  FracShootLeafElmAlloc2Litr(ielmp,k_woody_litr)=FracShootLeafElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPSHW

  FracRootStalkElmAlloc2Litr(ielmn,k_woody_litr)=FracRootStalkElmAlloc2Litr(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNRTW
  FracRootStalkElmAlloc2Litr(ielmp,k_woody_litr)=FracRootStalkElmAlloc2Litr(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPRTW

  FracRootElmAlloc2Litr(ielmn,k_woody_litr)=FracRootElmAlloc2Litr(ielmc,k_woody_litr)*RootrNC_pft(NZ)/CNRTW
  FracRootElmAlloc2Litr(ielmp,k_woody_litr)=FracRootElmAlloc2Litr(ielmc,k_woody_litr)*RootrPC_pft(NZ)/CPRTW

  DO NE=2,NumPlantChemElms
    FracShootStalkElmAlloc2Litr(NE,k_fine_litr)=1.0_r8-FracShootStalkElmAlloc2Litr(NE,k_woody_litr)
    FracShootLeafElmAlloc2Litr(NE,k_fine_litr)=1.0_r8-FracShootLeafElmAlloc2Litr(NE,k_woody_litr)
    FracRootStalkElmAlloc2Litr(NE,k_fine_litr)=1.0_r8-FracRootStalkElmAlloc2Litr(NE,k_woody_litr)
    FracRootElmAlloc2Litr(NE,k_fine_litr)=1.0_r8-FracRootElmAlloc2Litr(NE,k_woody_litr)
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
  TKCM=TKC_pft(NZ)+TempOffset_pft(NZ)

  TFN5=calc_plant_maint_tempf(TKCM)  
  D7: DO L=NU,MaxNumRootLays
    TKSM=TKS_vr(L)+TempOffset_pft(NZ)
    TFN6_vr(L)=calc_plant_maint_tempf(TKSM)  
  ENDDO D7
!
!     PRIMARY ROOT NUMBER
!
!     RootBiomCPerPlant_pft=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     RootPrimeAxsNum=multiplier for number of primary root axes
!
  RootBiomCPerPlant_pft(NZ)=AMAX1(0.999992087_r8*RootBiomCPerPlant_pft(NZ),&
    RootElms_pft(ielmc,NZ)/PlantPopulation_pft(NZ))
  RootPrimeAxsNum=AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**0.667_r8)*PlantPopulation_pft(NZ)
!
!     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
!     FROM CANOPY TURGOR
!
!     WFNS=turgor expansion,extension function
!     PSICanopyTurg_pft,PSIMin4OrganExtens=current,minimum canopy turgor potl for expansion,extension
!     Stomata_Stress=stomatal resistance function of canopy turgor
!     PSICanopy_pft=canopy water potential
!     WaterStress4Groth=growth function of canopy water potential
!     CanTurgPSIFun4Expans=expansion,extension function of canopy water potential
!
  WFNS=AMIN1(1.0_r8,AZMAX1(PSICanopyTurg_pft(NZ)-PSIMin4OrganExtens))

  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    !bryophyte, no turgor
    Stomata_Stress    = 1.0_r8
    WaterStress4Groth = EXP(0.05_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8))
  ELSE
    !others
    Stomata_Stress    = EXP(RCS_pft(NZ)*PSICanopyTurg_pft(NZ))
    WaterStress4Groth = EXP(0.10_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8))
  ENDIF
!  write(119,*)I*1000+J,Stomata_Stress,RCS_pft(NZ),PSICanopyTurg_pft(NZ),PSICanopy_pft(NZ),plt_ew%PSICanopyOsmo_pft(NZ)
  CanTurgPSIFun4Expans            = fRespWatSens(WFNS,iPlantRootProfile_pft(NZ))
  plt_biom%StomatalStress_pft(NZ) = Stomata_Stress
  end associate
  end subroutine StagePlantForGrowth
!------------------------------------------------------------------------------------------

  subroutine ComputeTotalBiom(I,J,NZ)

  integer, intent(in) :: I,J,NZ
  
  integer :: L,N
!     begin_execution
  associate(                                                   &
    Myco_pft                 => plt_morph%Myco_pft,                &  !input
    MaxSoiL4Root_pft       => plt_morph%MaxSoiL4Root_pft,      &  !input
    NU                     => plt_site%NU,                     &  !input
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr, &  !input
    PopuRootMycoC_pvr      => plt_biom% PopuRootMycoC_pvr      &  !inoput
  )

  call SumPlantBiom(I,J,NZ,'computotb')

  if(lverb)write(*,*)'add the nonstrucal components'
  D3451: DO N=1,Myco_pft(NZ)
    DO  L=NU,MaxSoiL4Root_pft(NZ)
      PopuRootMycoC_pvr(N,L,NZ)=PopuRootMycoC_pvr(N,L,NZ)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
    enddo
  ENDDO D3451

  end associate

  end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

  subroutine AccumulateStates(I,J,NZ,CanopyN2Fix_pft)
  
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), optional, intent(in) :: CanopyN2Fix_pft(JP1)
  integer :: L,NR,N,NE,NB

!     begin_execution
  associate(                                                         &
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch, &  !input
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch, &  !input
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,      &  !input
    CanopyStalkArea_lbrch     => plt_morph%CanopyStalkArea_lbrch,    &  !input
    EarStrutElms_brch         => plt_biom%EarStrutElms_brch,         &  !input
    GrainStrutElms_brch       => plt_biom%GrainStrutElms_brch,       &  !input
    HuskStrutElms_brch        => plt_biom%HuskStrutElms_brch,        &  !input
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch,    &  !input
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,         &  !input
    NU                        => plt_site%NU,                        &  !input
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,        &  !input
    RootH2PO4Uptake_pft       => plt_rbgc%RootH2PO4Uptake_pft,       &  !input
    RootHPO4Uptake_pft        => plt_rbgc%RootHPO4Uptake_pft,        &  !input
    RootMycoExudElms_pft      => plt_rbgc%RootMycoExudElms_pft,      &  !input
    RootN2Fix_pft             => plt_rbgc%RootN2Fix_pft,             &  !input
    RootNH4Uptake_pft         => plt_rbgc%RootNH4Uptake_pft,         &  !input
    RootNO3Uptake_pft         => plt_rbgc%RootNO3Uptake_pft,         &  !input
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr,   &  !input
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr,   &  !input
    SeedNumSet_brch           => plt_morph%SeedNumSet_brch,          &  !input
    ShootStrutElms_brch       => plt_biom%ShootStrutElms_brch,       &  !input
    StalkLiveBiomassC_brch    => plt_biom%StalkLiveBiomassC_brch,    &  !input
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch,        &  !input
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch,       &  !input
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,       &  !input
    CanopyNodulNonstElms_pft  => plt_biom%CanopyNodulNonstElms_pft,  &  !inoput
    CanopyStemAreaZ_pft       => plt_morph%CanopyStemAreaZ_pft,      &  !inoput
    LeafAreaLive_brch         => plt_morph%LeafAreaLive_brch,        &  !inoput
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch,        &  !inoput
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch,      &  !inoput
    PlantExudElm_CumYr_pft    => plt_rbgc%PlantExudElm_CumYr_pft,    &  !inoput
    PlantN2Fix_CumYr_pft      => plt_bgcr%PlantN2Fix_CumYr_pft,      &  !inoput
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft,  &  !inoput
    RootUptk_N_CumYr_pft      => plt_rbgc%RootUptk_N_CumYr_pft,      &  !inoput
    RootUptk_P_CumYr_pft      => plt_rbgc%RootUptk_P_CumYr_pft,      &  !inoput
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft,       &  !output
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,      &  !output
    CanopyNodulElms_pft       => plt_biom%CanopyNodulElms_pft,       &  !output
    CanopyNonstElms_pft       => plt_biom%CanopyNonstElms_pft,       &  !output
    CanopySeedNum_pft         => plt_morph%CanopySeedNum_pft,        &  !output
    CanopyStalkC_pft          => plt_biom%CanopyStalkC_pft,          &  !output
    CanopyStemArea_pft        => plt_morph%CanopyStemArea_pft,       &  !output
    EarStrutElms_pft          => plt_biom%EarStrutElms_pft,          &  !output
    GrainStrutElms_pft        => plt_biom%GrainStrutElms_pft,        &  !output
    HuskStrutElms_pft         => plt_biom%HuskStrutElms_pft,         &  !output
    LeafStrutElms_pft         => plt_biom%LeafStrutElms_pft,         &  !output
    PetoleStrutElms_pft       => plt_biom%PetoleStrutElms_pft,       &  !output
    RootNodulElms_pft         => plt_biom%RootNodulElms_pft,         &  !output
    ShootStrutElms_pft        => plt_biom%ShootStrutElms_pft,        &  !output
    StalkRsrvElms_pft         => plt_biom%StalkRsrvElms_pft,         &  !output
    StalkStrutElms_pft        => plt_biom%StalkStrutElms_pft         &  !output
  )
!
!     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
!
!
  call SumPlantBiom(I,J,NZ,'Accumstates')

  DO NB=1,NumOfBranches_pft(NZ)    
    if(isclose(LeafPetolBiomassC_brch(NB,NZ),0._r8))then
      LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
      PetoleStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
      LeafAreaLive_brch(NB,NZ)                     = 0._r8
    endif
 ENDDO  
  DO NE=1,NumPlantChemElms
    CanopyNonstElms_pft(NE,NZ) = AZMAX1(sum(CanopyNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    StalkRsrvElms_pft(NE,NZ)   = AZMAX1(sum(StalkRsrvElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    ShootStrutElms_pft(NE,NZ)  = AZMAX1(sum(ShootStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    PetoleStrutElms_pft(NE,NZ) = AZMAX1(sum(PetoleStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    StalkStrutElms_pft(NE,NZ)  = AZMAX1(sum(StalkStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    LeafStrutElms_pft(NE,NZ)   = AZMAX1(sum(LeafStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    HuskStrutElms_pft(NE,NZ)   = AZMAX1(sum(HuskStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    GrainStrutElms_pft(NE,NZ)  = AZMAX1(sum(GrainStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    EarStrutElms_pft(NE,NZ)    = AZMAX1(sum(EarStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    !root state variables
    !sum structural biomass
  ENDDO

  CanopyStalkC_pft(NZ)     = AZMAX1(sum(StalkLiveBiomassC_brch(1:NumOfBranches_pft(NZ),NZ)))
  CanopyLeafShethC_pft(NZ) = AZMAX1(sum(LeafPetolBiomassC_brch(1:NumOfBranches_pft(NZ),NZ)))
  CanopySeedNum_pft(NZ)    = AZMAX1(sum(SeedNumSet_brch(1:NumOfBranches_pft(NZ),NZ)))
  CanopyLeafArea_pft(NZ)   = AZMAX1(sum(LeafAreaLive_brch(1:NumOfBranches_pft(NZ),NZ)))
  CanopyStemArea_pft(NZ)   = AZMAX1(sum(CanopyStalkArea_lbrch(1:NumOfCanopyLayers1,1:NumOfBranches_pft(NZ),NZ)))

  
  DO NB=1,NumOfBranches_pft(NZ)        
    DO L=1,NumOfCanopyLayers1
      CanopyStemAreaZ_pft(L,NZ)=CanopyStemAreaZ_pft(L,NZ)+CanopyStalkArea_lbrch(L,NB,NZ)
    ENDDO
  ENDDO

!
!     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
!
!     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
  IF(is_plant_N2fix(iPlantNfixType_pft(NZ)))THEN
    IF(is_canopy_N2fix(iPlantNfixType_pft(NZ)))THEN
      DO NE=1,NumPlantChemElms
        D7950: DO NB=1,NumOfBranches_pft(NZ)
          CanopyNodulNonstElms_pft(NE,NZ)=CanopyNodulNonstElms_pft(NE,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ)
        ENDDO D7950
        CanopyNodulElms_pft(NE,NZ)=sum(CanopyNodulStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))+ &
          sum(CanopyNodulNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
      ENDDO
    ELSEIF(is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
      DO NE=1,NumPlantChemElms
        RootNodulElms_pft(NE,NZ)=sum(RootNodulStrutElms_rpvr(NE,NU:MaxSoiL4Root_pft(NZ),NZ))+&
          sum(RootNodulNonstElms_rpvr(NE,NU:MaxSoiL4Root_pft(NZ),NZ))
      ENDDO
    ENDIF
  ENDIF
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

  PlantN2Fix_CumYr_pft(NZ)=PlantN2Fix_CumYr_pft(NZ)+RootN2Fix_pft(NZ)
  if(present(CanopyN2Fix_pft))PlantN2Fix_CumYr_pft(NZ)=PlantN2Fix_CumYr_pft(NZ)+CanopyN2Fix_pft(NZ)
  end associate
  end subroutine AccumulateStates

!------------------------------------------------------------------------------------------

  subroutine ZeroPlantStates(I,J,NZ)  
  !
  !for annual plants
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NE
  DO NE=1,NumPlantChemElms
    plt_biom%CanopyNonstElms_pft(NE,NZ)=0._r8
    plt_biom%StalkRsrvElms_pft(NE,NZ)=0._r8
    plt_biom%ShootStrutElms_pft(NE,NZ)=0._r8
    plt_biom%PetoleStrutElms_pft(NE,NZ)=0._r8
    plt_biom%StalkStrutElms_pft(NE,NZ)=0._r8
    plt_biom%LeafStrutElms_pft(NE,NZ)=0._r8
    plt_biom%HuskStrutElms_pft(NE,NZ)=0._r8
    plt_biom%GrainStrutElms_pft(NE,NZ)=0._r8
    plt_biom%EarStrutElms_pft(NE,NZ)=0._r8
    plt_biom%CanopyNodulElms_pft(NE,NZ)=0._r8
    plt_biom%RootNodulElms_pft(NE,NZ)=0._r8
    plt_rbgc%PlantRootSoilElmNetX_pft(NE,NZ)=0._r8
  ENDDO

  plt_biom%CanopyStalkC_pft(NZ)=0._r8
  plt_biom%CanopyLeafShethC_pft(NZ) =0._r8
  plt_morph%CanopySeedNum_pft(NZ) =0._r8
  plt_morph%CanopyLeafArea_pft(NZ)=0._r8
  plt_morph%CanopyStemArea_pft(NZ)=0._r8
  plt_morph%CanopyStemAreaZ_pft(1:NumOfCanopyLayers1,NZ)=0._r8

  end subroutine ZeroPlantStates

end module grosubsMod
